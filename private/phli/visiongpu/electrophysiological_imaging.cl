/* Defined by preprocessor options:
 *   NLPOINTS
 *   NRPOINTS
 *   BUFFERSIZE
 *   NUMELECTRODES
 *   MAXSPIKES
 *   LSIZE
 */

__kernel void electrophysiological_imaging(
	__global const short8* samples, 
    __constant const short* cellnums, __constant const short* spikenums, __constant const short* spikeinds, 
    __global const short* spikelist,
	__global float* baseline, __global short8* waveform, __global short8* error,
	int sampleindex,
	__local short8* lsamples)
{
	int groupid = get_group_id(0);
	if (groupid >= NUMELECTRODES/8) return;
	
	int numgroups0 = get_num_groups(0);
	int lid = get_local_id(0);


	// Collaboratively load samples data for this group's electrodes into local memory	
	int samplenum = lid;
	for (int i = 0; i < BUFFERSIZE/LSIZE; i++) {
		lsamples[samplenum] = samples[samplenum*NUMELECTRODES/8 + groupid];
		samplenum += LSIZE;
	}
	if (samplenum < BUFFERSIZE) lsamples[samplenum] = samples[samplenum*NUMELECTRODES/8 + groupid];
	
	// Wait for local memory to be filled with sample data
	barrier(CLK_LOCAL_MEM_FENCE);


	// Calculate cellindex
	int gid1 = get_global_id(1);
	int lsize = get_local_size(0);
	int cellindex = lid + gid1*lsize;

	// This seems to run fine whether before or after the barrier on creampuff with the test data.  
	// Not clear if it will cause crashes if put before the barrier on some systems due to some threads 
	// in a warp/wavefront not making it to the barrier while others do.
	if (cellindex > NUMNEURONS) return; 

	// Get spike info
	short numspikes = spikenums[cellindex];
	short spike_start_index = spikeinds[cellindex];
	
	// Get cellnum, calculate waveform index
	short cellnum = cellnums[cellindex];
	int waveform_cell_index = cellnum * (NLPOINTS+1+NRPOINTS) * NUMELECTRODES/8;
	
	
	// Calculate EIs!
	short8 point8;
	for (int p = -NLPOINTS; p <= NRPOINTS; p++) {
		point8 = 0;
		for (int s = 0; s < numspikes; s++) {
			short spike = spikelist[spike_start_index+s];
			point8 += lsamples[spike+p];
		}
		
		// This line is killing it; out of resources OPENCL says...
		waveform[waveform_cell_index + (p+NLPOINTS)*NUMELECTRODES/8 + groupid] = point8;
	}

	
	
	
	
	// Hacked sanity checks; just use waveform/error as dumping ground.  
	// Easier to read if waveform/error are just short* not short8*
//	if (groupid != 0 || lid != 0 || gid1 != 0) return;
//	waveform[sampleindex / BUFFERSIZE] = numspikes;	
//	error[sampleindex / BUFFERSIZE] = cellnum;
//	waveform[sampleindex / BUFFERSIZE] = lsamples[0].s2;
}