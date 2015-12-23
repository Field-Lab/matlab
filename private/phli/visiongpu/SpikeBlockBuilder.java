package edu.ucsc.neurobiology.vision.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import edu.ucsc.neurobiology.vision.io.NeuronFile;

/**
 * 
 * @author Peter H. Li, The Salk Institute, peterli@salk.edu
 */
public class SpikeBlockBuilder {
	private final NeuronFile neuronFile;
	private final int nlpoints, nrpoints;
	private final int blockSize, nBlocks, spikesPerBlock;
	private final int[] ids;
	private final ArrayList<LinkedList<Integer>> spikeTimes;
	private int sampleIndex = 0;
	
	public SpikeBlockBuilder(NeuronFile nf, int spikesToUse, int nlpoints, int nrpoints, int blockSize) throws IOException {
		this.neuronFile = nf;
		this.nlpoints = nlpoints;
		this.nrpoints = nrpoints;
		
		this.blockSize = blockSize;
		this.nBlocks = 1 + (nf.getNumberOfSamples() - 1) / blockSize;

		int spikesPerBlock = Math.round((float) spikesToUse / nBlocks);
		if (spikesToUse < 1) spikesPerBlock = Integer.MAX_VALUE;
		if (spikesPerBlock < 1) spikesPerBlock = 1;
		this.spikesPerBlock = spikesPerBlock;
		
		this.ids = neuronFile.getIDList();
		this.spikeTimes = new ArrayList<LinkedList<Integer>>(neuronFile.getNumberOfNeurons());
		for (int id : ids) {
			int[] times = neuronFile.getSpikeTimes(id);
			LinkedList<Integer> spikeList = new LinkedList<Integer>();
			for (int t : times) spikeList.add(t);
			spikeTimes.add(spikeList);
		}
	}
	
	public int getBlockSize() {	return blockSize; }
	public int getNumNeurons() { return neuronFile.getNumberOfNeurons(); }
	
	@Override
	public String toString() {
		return "SpikeBlockBuilder:\n\n" +
				"neuronFile:\n" + neuronFile + "\n" +
				"nlpoints, nrpoints: " + nlpoints + "," + nrpoints + "\n" +
				"blockSize: " + blockSize + "\n" +
				"nBlocks: " + nBlocks + "\n" +
				"spikesPerBlock: " + spikesPerBlock;
	}
	
	public RawSpikeBlock nextRaw() {
		LinkedList<NeuronSpikeBlock> block = new LinkedList<NeuronSpikeBlock>();
		
		for (int n = 0; n < ids.length; n++) {
			LinkedList<Integer> spikeList = spikeTimes.get(n);
			LinkedList<Integer> inBlock = new LinkedList<Integer>();
			int time;

			// Skip past spikes that are too close to beginning of block
			while (!spikeList.isEmpty()) {
				time = spikeList.removeFirst();
				if (time >= sampleIndex + nlpoints) {
					spikeList.addFirst(time);
					break; 
				}
			}

			// Collect spikes within block
			while (!spikeList.isEmpty()) {
				time = spikeList.removeFirst();
				if (time >= sampleIndex + blockSize - nrpoints) {
					spikeList.addFirst(time);
					break;
				}
				inBlock.add(time);
			}
			
			// Skip some spikes to get the right number
			int nInBlock = inBlock.size();
			int step = Math.round((float) nInBlock / spikesPerBlock);
			if (step < 1) step = 1;
			
			// Collect selected spikes
			LinkedList<Short> selected = new LinkedList<Short>();
			for (int s = 0; s < nInBlock; s += step) selected.add((short) (inBlock.get(s) - sampleIndex));
			
			block.add(new NeuronSpikeBlock(ids[n], n, selected));
		}

		// Sort the NeuronSpikeBlocks in order of number of spikes (important for efficient branching on GPU)
		Collections.sort(block);

		RawSpikeBlock rawBlock = new RawSpikeBlock(block);
		sampleIndex += blockSize;
		return rawBlock;
	}

	
	/**
	 * Collects NeuronSpikeBlocks for all neurons
	 */
	class RawSpikeBlock {
		private final List<NeuronSpikeBlock> neuronSpikeBlocks;
		public final int nNeurons;
		
		private HashMap<Integer,Integer> idToIndexHash = null;
		private int numSpikes = -1;
		
		public RawSpikeBlock(List<NeuronSpikeBlock> blocks) {
			this.neuronSpikeBlocks = blocks;
			this.nNeurons = neuronSpikeBlocks.size();
		}
		
		private HashMap<Integer,Integer>getIDToIndexHash() {
			if (idToIndexHash == null) initHash();
			return idToIndexHash;
		}
		
		private void initHash() {
			idToIndexHash = new HashMap<Integer,Integer>(nNeurons);
			for (int i = 0; i < nNeurons; i++)
				idToIndexHash.put(neuronSpikeBlocks.get(i).cellID, i);
		}
		
		public NeuronSpikeBlock getID(int id) { return get(getIDToIndexHash().get(id)); }
		public NeuronSpikeBlock get(int index) { return neuronSpikeBlocks.get(index); }
		
		public short[] getCellNums() {
			short[] cellNums = new short[nNeurons];
			for (int i = 0; i < nNeurons; i++) cellNums[i] = (short) neuronSpikeBlocks.get(i).cellNum;
			return cellNums;
		}
		
		public short[] getSpikeNums() {
			short[] spikeNums = new short[nNeurons];
			for (int i = 0; i < nNeurons; i++) spikeNums[i] = neuronSpikeBlocks.get(i).numSpikes;
			return spikeNums;
		}
		
		public short[] getSpikeIndices() {
			short[] spikeIndices = new short[nNeurons];
			for (int i = 1; i < nNeurons; i++) spikeIndices[i] = (short) (neuronSpikeBlocks.get(i-1).numSpikes + spikeIndices[i-1]);
			return spikeIndices;
		}
		
		public int getNumSpikes() {
			if (numSpikes >= 0) return numSpikes;

			numSpikes = 0;
			for (NeuronSpikeBlock nsb : neuronSpikeBlocks) numSpikes += nsb.numSpikes;
			return numSpikes;
		}
		
		public short[] getSpikeList() {
			short[] spikeList = new short[getNumSpikes()];
			int i = 0;
			for (NeuronSpikeBlock nsb : neuronSpikeBlocks) {
				for (Short spike : nsb.spikeList) 
					spikeList[i++] = spike;
			}
			return spikeList;
		}
		
	}
		
	
	/**
	 * Holds the list of spikes that fall within a given block for one neuron.
	 */
	class NeuronSpikeBlock implements Comparable<NeuronSpikeBlock> {
		public final int cellID;
		public final int cellNum;
		public final LinkedList<Short> spikeList;
		public final short numSpikes;
		
		public NeuronSpikeBlock(int cellID, int cellNum, LinkedList<Short> selected) {
			this.cellID = cellID;
			this.cellNum = cellNum;
			this.spikeList = selected;
			this.numSpikes = (short) selected.size();
		}

		public int compareTo(NeuronSpikeBlock other) {
			return other.numSpikes - this.numSpikes;
		}

	}
	
	
}