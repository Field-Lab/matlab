package edu.ucsc.neurobiology.vision.analysis;

import java.io.IOException;
import java.util.Arrays;

import org.bridj.Pointer;

import com.nativelibs4java.opencl.*;
import com.nativelibs4java.opencl.CLKernel.LocalSize;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.util.IOUtils;

import edu.ucsc.neurobiology.vision.analysis.SpikeBlockQueueCL.SpikeBlockCL;
import edu.ucsc.neurobiology.vision.io.AsyncBufferCL;
import edu.ucsc.neurobiology.vision.io.NeuronFile;
import edu.ucsc.neurobiology.vision.io.RawDataHeader;

public class EICalculatorCL {
	private static int CSHORT_SIZE = 2; // Couldn't find this in Bridj but it must be there...
//	private static final int MAX_WORKITEMS_PER_GROUP = 256; // ATI Evergreen series (Radeon HD 5xxx, 6xxx)
	private static final int MAX_WORKITEMS_PER_GROUP = 512; // NVIDIA can go to 1024?  But 512 maybe more efficient
	private static final int PREFERRED_ELECTRODES_PER_GROUP = 8; // 8 shorts fits well into 128-bit memory bus
	
	private final int nElectrodes, nNeurons;
	private final int nPoints;
	
	private final int nNeuronBlocks;
	private final int nGroups, nWorkitems0;
	
	private final CLQueue clQueue;
	private final CLKernel eiKernel;
	private final CLBuffer<Float> movingBaseline;
	private final CLBuffer<Short> waveformAccumulator, waveformErrorAccumulator;
	private final LocalSize localSize;

	
	public EICalculatorCL(NeuronFile neuronFile, RawDataHeader rawDataHeader, 
			int nlpoints, int nrpoints, 
			int bufferSize, int maxSpikesPerCellPerBlock, 
			CLContext context, CLQueue clQueue) throws IOException {
		
		this.nNeurons = neuronFile.getNumberOfNeurons();
		this.nElectrodes = rawDataHeader.getNumberOfElectrodes();
		this.nPoints = nlpoints + 1 + nrpoints;

		this.nNeuronBlocks = ((nNeurons - 1) / MAX_WORKITEMS_PER_GROUP) + 1;
		this.nGroups = ((nElectrodes - 1) / PREFERRED_ELECTRODES_PER_GROUP) + 1;
		this.nWorkitems0 = nGroups * MAX_WORKITEMS_PER_GROUP;
		this.localSize = new LocalSize(CSHORT_SIZE * PREFERRED_ELECTRODES_PER_GROUP * bufferSize);
		
		this.clQueue = clQueue;
		String src = IOUtils.readText(EICalculatorCL.class.getResource("electrophysiological_imaging.cl"));
		CLProgram program = context.createProgram(src);	

		program.addBuildOption("-DNLPOINTS=" + nlpoints);
		program.addBuildOption("-DNRPOINTS=" + nrpoints);
		program.addBuildOption("-DBUFFERSIZE=" + bufferSize);
		program.addBuildOption("-DNUMELECTRODES=" + nElectrodes);
		program.addBuildOption("-DNUMNEURONS=" + nNeurons);
		program.addBuildOption("-DMAXSPIKES=" + maxSpikesPerCellPerBlock);
		program.addBuildOption("-DLSIZE=" + MAX_WORKITEMS_PER_GROUP);
//		program.addBuildOption("-cl-fast-relaxed-math");
		eiKernel = program.createKernel("electrophysiological_imaging");
		
		// Create GPU buffers for EI calculation
        movingBaseline           = context.createBuffer(Usage.InputOutput, Float.class, nElectrodes);
        waveformAccumulator      = context.createBuffer(Usage.InputOutput, Short.class, nElectrodes * nPoints * neuronFile.getNumberOfNeurons());
        waveformErrorAccumulator = context.createBuffer(Usage.InputOutput, Short.class, nElectrodes * nPoints * neuronFile.getNumberOfNeurons());
	}
	
	
	public void processSamples(int sampleIndex, AsyncBufferCL<Short> sampleBuffer, SpikeBlockCL spikeBlock) {
		eiKernel.setArgs(
				sampleBuffer.getBuffer(), 
				spikeBlock.cellNums.getBuffer(), spikeBlock.spikeNums.getBuffer(), spikeBlock.spikeIndices.getBuffer(), 
				spikeBlock.spikeList.getBuffer(),
				movingBaseline, waveformAccumulator, waveformErrorAccumulator, 
				sampleIndex,
				localSize);

		// FIXME: Having some weird issues with spikeBlock buffers not being ready to read from when reached in the kernel.
		// Have to look into more.
		sampleBuffer.waitFor();
		spikeBlock.waitFor();
		
		CLEvent eiEvt = eiKernel.enqueueNDRange(clQueue, new int[]{ nWorkitems0, nNeuronBlocks }, new int[]{ MAX_WORKITEMS_PER_GROUP, 1 });
		sampleBuffer.withEvt(eiEvt);
		spikeBlock.withEvt(eiEvt);
	}
	
	
	public void finishCalculation() {
		// Sanity check that we got to the right spike number
//		Pointer<Integer> spikeIndexPtr = spikeIndex.read(clQueue);
//		int spikeIndexFinal = spikeIndexPtr.getInt();
//		System.out.println("Final spike index: " + spikeIndexFinal);
//		System.out.println("Final spike time: " + (double) allSpikes[0][spikeIndexFinal] / header.getSamplingFrequency() + " s");
//		System.out.println("Penultimate spike time: " + (double) allSpikes[0][spikeIndexFinal-1] / header.getSamplingFrequency() + " s");
		
		Pointer<Short> waveformPtr      = waveformAccumulator.read(clQueue);
		Pointer<Short> waveformErrorPtr = waveformErrorAccumulator.read(clQueue);	
//		Pointer<Float> baselinePtr = movingBaseline.read(clQueue);

		short[] hacks = waveformPtr.getShorts(100);
		System.out.println(Arrays.toString(hacks));
		
		hacks = waveformErrorPtr.getShorts(100);
		System.out.println(Arrays.toString(hacks));		
	}
	
}
