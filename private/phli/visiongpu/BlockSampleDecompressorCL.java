package edu.ucsc.neurobiology.vision.io;

import java.io.File;
import java.io.FileInputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.ListIterator;

import org.bridj.Pointer;

import com.nativelibs4java.opencl.*;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.util.IOUtils;

import edu.ucsc.neurobiology.vision.util.VisionParams;

/**
 * Wraps a raw data file and the OpenCL kernel for decompressing raw data on the GPU.
 * 
 * Need to make this work for directories of raw data files as well.
 * 
 * @author Peter H. Li, The Salk Institute
 * @owner Peter H. Li, The Salk Institute
 */
public class BlockSampleDecompressorCL {
	private static final int BYTES_IN_PER_OP = 3;
	private static final int SHORTS_OUT_PER_OP = 2;
	
	private final CLContext context;
	private final CLQueue queue;
	private final CLKernel decompressSamplesKernel;
	private final RawDataHeader header;
	private InputStream[] files;
	private int file = 0;
	
	public static BlockSampleDecompressorCL create(CLContext context, CLQueue queue, String rawDataFilename) throws IOException {
		return create(context, queue, new File(rawDataFilename));
	}
	
	public static BlockSampleDecompressorCL create(CLContext context, CLQueue queue, File rawDataFile) throws IOException {
		InputStream[] rawDataStreams = null;
		if (rawDataFile.isFile()) {
			FileInputStream rawDataStream = new FileInputStream(rawDataFile);
			rawDataStreams = new InputStream[]{ rawDataStream };
		} else if(rawDataFile.isDirectory()) {
			File[] rawDataFiles = rawDataFile.listFiles(new FilenameFilter() {
				public boolean accept(File dir, String name) {
					return name.endsWith(VisionParams.BIN_FILE_EXTENSION_512);
				}
			});
			List<File> rawDataFileList = Arrays.asList(rawDataFiles);
			Collections.sort(rawDataFileList);
			rawDataStreams = new InputStream[rawDataFileList.size()];
			for (int i = 0; i < rawDataStreams.length; i++) rawDataStreams[i] = new FileInputStream(rawDataFileList.get(i));
		}
		return new BlockSampleDecompressorCL(context, queue, rawDataStreams);
	}	
	
	public BlockSampleDecompressorCL(CLContext context, CLQueue queue, InputStream[] rawDataStreams) throws IOException {
		this.context = context;		
		this.queue = queue;
		
		this.files = rawDataStreams;
		this.header = RawDataHeader.extractHeader(files[0]);
		
		// Read the program sources and compile them into kernel
		String src = IOUtils.readText(BlockSampleDecompressorCL.class.getResource("decompress_samples.cl"));
		CLProgram program = context.createProgram(src);	
		if (header.getNumberOfElectrodes() % 2 == 0)
			this.decompressSamplesKernel = program.createKernel("decompress_samples");
		else 
			this.decompressSamplesKernel = program.createKernel("decompress_samples_odd");
	}

	
	public RawDataHeader getHeader() { return header; }
	

	public CLBuffer<Short> createOutBuffer(int samples) {
		return context.createBuffer(Usage.InputOutput, Short.class, samples * header.getNumberOfElectrodes());
	}
	
	public CLBuffer<Byte> createInBuffer(int samples) {
		return context.createBuffer(Usage.Input, Byte.class, samples * header.getSampleSize());
	}
	
	
	public AsyncBufferCL<Short> decompress(int samples) throws IOException {
		CLBuffer<Byte> in = createInBuffer(samples);
		CLBuffer<Short> out = createOutBuffer(samples);
		CLEvent decompressEvt = decompress(in, out);
		
		if (decompressEvt == null) return null;
		in.release();
		return new AsyncBufferCL<Short>(out).withEvt(decompressEvt);
	}
	
	
	/**
	 * The in and out buffers must be appropriate sizes otherwise behavior is undefined.
	 * 
	 * @param in
	 * @param out
	 * @return
	 * @throws IOException
	 */
	public CLEvent decompress(CLBuffer<Byte> in, CLBuffer<Short> out) throws IOException {
		int nBytes = (int) in.getByteCount();
		byte[] b = new byte[nBytes];
		
		int totRead = 0;
		while (totRead < nBytes) {
			// We are all out of data
			if (file == files.length) {
				if (totRead > 0) {
					Arrays.fill(b, totRead, nBytes, (byte) 0);
					break;
				} else {
					return null;
				}
			}

			int thisRead = files[file].read(b, totRead, nBytes-totRead);

			// Out of data in current file, try the next one
			if (thisRead == -1) {
				file++;
				continue;
			}
			
			totRead += thisRead;
		}
		
		// Wrap bytes in ByteBuffer objects that JavaCL likes
		ByteBuffer bb = ByteBuffer.wrap(b);

		// Write bytes into GPU memory buffers
		CLEvent inEvt = in.writeBytes(getQueue(), 0, nBytes, bb, false);
		
		int nElectrodes = header.getNumberOfElectrodes();
		int ops;
		if (nElectrodes % 2 == 0) {
			ops = (int) nBytes / BYTES_IN_PER_OP;
			decompressSamplesKernel.setArgs(in, out, ops);
		} else {
			int samples = nBytes / header.getSampleSize();
			int ops_per_sample = (nElectrodes-1)/2 + 1;
			ops = samples * ops_per_sample;
			decompressSamplesKernel.setArgs(in, out, nElectrodes, ops);
		}
		return decompressSamplesKernel.enqueueNDRange(getQueue(), new int[] { ops }, inEvt);	
	}
	
	public CLEvent decompress(AsyncBufferCL<Byte> in, AsyncBufferCL<Short> out) throws IOException {
		in.waitFor();
		out.waitFor();
		return decompress(in.getBuffer(), out.getBuffer());
	}	
	
	
	public void close() throws IOException {
		for (InputStream file : files) file.close();
	}

	
	public CLQueue getQueue() {
		return queue;
	}

}
