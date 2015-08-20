function data = getData(request)
%bright bars moving right
dataBBMR= [2007-8-24-4	10	2	4	8	48.2458	0.4746	48.2118	0.3992	48.2525	0.8105	0	0
2007-8-24-4	11	1	4	16	48.1782	0.6345	48.2217	0.3242	48.6505	1.2637	0	0
2007-8-24-4	2	9	8	8	96.2906	0.7759	97.0162	1.1235	96.0080	1.4017	0	0
2007-8-24-4	6	3	8	8	96.3778	0.8529	96.6834	1.6250	95.8793	1.7374	0	0
2007-03-27-1	12	14	8	8	94.2439	1.4895	92.6371	1.6115	95.7292	1.5379	103.0478	26.2789
2007-03-27-1	15	3	8	8	93.3637	1.6917	92.8020	1.8534	95.1667	1.6961	91.1863	12.34785
2007-8-24-4	7	3	8	16	96.2764	0.8370	96.8241	1.0359	96.0624	1.9367	0	0
2007-03-27-1	18	1	8	16	93.4745	1.4377	90.9224	1.4474	93.5113	1.7181	91.2134	9.8567
2007-8-24-4	4	1	16	8	192.5311	5.2138	205.3461	19.8409	190.0693	16.9549	0	0
2007-03-27-1	16	4	16	8	184.6126	10.2005	176.4689	16.2248	190.8418	24.8675	194.7933	101.3310
2007-8-24-4	5	4	16	16	192.3259	3.5522	195.8780	5.8959	188.3950	7.8328	0	0
2007-8-24-4	3	23	16	16	191.3860	2.8556	195.3012	5.5379	189.7446	6.1453	0	0
2007-03-27-1	13	1	16	16	187.5700	4.4339	179.5899	6.3662	186.8455	6.7636	204.2783	82.6860
2007-03-27-1	19	1	16	16	183.9616	4.7369	177.5880	9.0639	187.5362	12.5737	199.0751	75.2832];
%dark bars moving right
dataDBMR = [2007-03-27-1	  12	2	8	8	93.1327	2.5905	93.9007	0.7951	94.3198	1.3409	96.5501	19.6723
2007-03-27-1	  15	4	8	8	90.8243	3.2959	92.7500	0.7907	93.8439	2.1854	90.4033	27.6925
2007-03-27-1	  18	4	8	16	91.0135	2.6286	92.2978	0.6092	91.9926	2.5153	93.3387	8.9170
2007-03-27-1	  16 	1	16	8	184.1164	24.5050	182.9856	5.3567	192.1703	24.8277	196.1818	50.5889
2007-03-27-1	  13	2	16	16	181.6240	7.7479	183.4271	3.5344	183.6619	10.6257	188.1424	42.8660
2007-03-27-1	  19	4	16	16	179.9152	10.0789	181.9082	5.0169	181.3347	13.4602	188.6278	33.3646
2007-8-24-4	  10	4	4	8	48.0354	1.2482	48.3788	0.1816	48.9651	3.5017	0	0
2007-8-24-4	  11	4	4	16	48.4958	0.4048	48.2571	0.2206	49.5562	4.1184	0	0
2007-8-24-4	  2	7	8	8	96.9170	1.0014	96.6715	0.5256	96.1238	1.7942	0	0
2007-8-24-4	  6	4	8	8	97.3693	1.4487	96.7375	0.4937	95.9965	1.9912	0	0
2007-8-24-4	  7	1	8	16	96.4753	0.7539	97.1396	0.4863	96.6462	2.4552	0	0
2007-8-24-4	  4	3	16	8	193.8624	8.0861	194.9085	4.2517	193.9645	24.2768	0	0
2007-8-24-4	  5 	2	16	16	193.5828	5.0253	193.7174	2.3587	193.6413	8.5866	0	0
2007-8-24-4	  3	22	16	16	194.4410	4.6275	193.8906	2.1418	198.3383	8.8379	0	0];

%bright bars moving left
dataBBML = [2007-03-27-1	  12	16	8	8	94.8607	1.4684	94.7311	1.9928	94.9947	1.5695	95.3229	13.4943
2007-03-27-1	  15	2	8	8	94.0237	1.8506	93.6187	2.3835	93.9831	1.6944	92.2614	11.1204
2007-03-27-1	  18	2	8	16	90.2729	2.7171	93.0857	0.8032	92.4858	2.2720	94.9005	6.0652
2007-03-27-1	  16 	3	16	8	187.1908	10.1278	190.6449	23.5806	187.3632	22.7869	187.8501	53.6012
2007-03-27-1	  13	10	16	16	185.6508	6.4955	186.5406	12.8605	184.8077	9.0203	208.9805	62.4933
2007-03-27-1	  19	2	16	16	184.1038	6.9355	188.7219	10.3075	186.3452	14.9657	194.1702	66.8869
2007-8-24-4	  10	3	4	8	48.0093	0.5485	48.1834	0.3826	48.3782	0.6291	0	0
2007-8-24-4	  11	2	4	16	47.9242	0.6165	48.4239	0.4343	48.6302	1.5909	0	0
2007-8-24-4	  2	19	8	8	95.9987	0.7318	95.3258	1.3285	96.8363	1.4306	0	0
2007-8-24-4	  6	1	8	8	95.9071	0.9401	95.7604	1.4990	97.0385	2.2269	0	0
2007-8-24-4	  7	4	8	16	95.6013	1.0300	95.8385	0.9679	96.7641	1.8774	0	0
2007-8-24-4	  4	2	16	8	191.5647	5.3579	192.1096	23.3616	188.3977	23.9724	0	0
2007-8-24-4	  5 	3	16	16	190.4741	3.9648	189.8284	6.7575	194.3592	8.8499	0	0
2007-8-24-4	  3	17	16	16	191.0115	3.4113	188.5061	6.3803	194.4348	8.8903	0	0];



dataDBML = [2007-03-27-1	  12	5	8	8	93.7029	3.3850	94.8136	0.6999	95.1391	1.9639	98.1572	15.6852
2007-03-27-1	  15	1	8	8	93.0909	3.8252	94.1667	0.8389	94.0491	2.044	97.7999	12.9781
2007-03-27-1	  18	3	8	16	90.2729	2.7171	93.0857	0.8032	92.4858	2.2720	94.9005	6.0652
2007-03-27-1	  16 	2	16	8	188.4690	20.4162	189.2267	8.4526	187.2130	22.7381	198.6514	51.6750
2007-03-27-1	  13	7	16	16	187.6019	11.2079	187.5061	4.7104	188.3484	10.7142	191.2378	36.6428
2007-03-27-1	  19	3	16	16	185.3288	13.9402	187.1187	4.5497	190.0643	16.7014	173.6851	32.7498
2007-8-24-4	  10	1	4	8	48.1388	0.5865	48.1169	0.1886	48.2420	0.7454	0	0
2007-8-24-4	  11	3	4	16	47.6074	0.4337	48.3481	0.2361	47.9567	1.3286	0	0
2007-8-24-4	  2	17	8	8	95.2699	1.0294	94.6568	1.2803	96.3501	1.9308	0	0
2007-8-24-4	  6	2	8	8	95.4485	1.1108	94.6802	1.0343	96.0173	1.6618	0	0
2007-8-24-4	  7	2	8	16	95.2436	0.9315	96.0347	0.5052	95.9390	2.0378	0	0
2007-8-24-4	  4	4	16	8	189.0281	6.3750	188.7123	3.9701	199.5215	16.6594	0	0
2007-8-24-4	  5 	1	16	16	188.8079	5.5961	190.4295	8.3233	187.4636	9.6448	0	0
2007-8-24-4	  3	2	16	16	189.1698	5.1273	188.8348	8.1405	189.6262	8.7575	0	0];


dataTau = [2007-03-27-1	12	14	8	8	98.0719	16.4754	95.2576	14.6778	98.4861	40.9538	91.0562	17.6953	0.1
2007-03-27-1	15	3	8	8	97.9582	14.76	96.7209	14.7166	90.6520	17.8048	88.0431	15.9372	0.1
2007-03-27-1	18	1	8	16	97.0594	9.6867	94.3343	9.4708	118.3704	122.9704	83.3876	12.1662	0.1
2007-03-27-1	16	4	16	8	195.7579	41.4464	210.9673	73.0297	188.7395	46.4675	197.7890	49.2968	0.1
2007-03-27-1	13	1	16	16	188.8643	41.5041	228.2198	41.6297	190.6702	43.2921	208.5265	86.0240	0.1
2007-03-27-1	19	1	16	16	196.2794	44.2047	236.4097	77.1543	198.4872	48.7742	192.6108	36.1036	0.1
2007-03-27-1	12	14	8	8	96.6938	9.9301	93.3894	8.3042	93.8442	15.0943	93.8489	17.8692	0.001
2007-03-27-1	15	3	8	8	96.0864	11.6463	92.7042	10.8765	96.5552	15.6884	94.1471	15.3331	0.001
2007-03-27-1	18	1	8	16	94.9634	13.5026	90.6465	10.5620	99.2145	14.1204	97.3337	14.9585	0.001
2007-03-27-1	16	4	16	8	187.3365	26.8035	190.4633	29.6493	192.4071	31.6995	195.0860	31.5828	0.001
2007-03-27-1	13	1	16	16	194.7054	20.4767	185.5863	22.2662	192.5745	25.9934	182.1464	32.8034	0.001
2007-03-27-1	19	1	16	16	188.9796	20.0535	188.8244	25.703	181.0396	25.3133	192.9131	33.33721	0.001
2007-03-27-1	12	14	8	8	94.2439	1.4895	92.6371	1.6115	95.7292	1.5379	103.0478	26.2789	0.01
2007-03-27-1	15	3	8	8	93.3637	1.6917	92.8020	1.8534	95.1667	1.6961	91.1863	12.34785	0.01
2007-03-27-1	18	1	8	16	93.4745	1.4377	90.9224	1.4474	93.5113	1.7181	91.2134	9.8567	0.01
2007-03-27-1	16	4	16	8	184.6126	10.2005	176.4689	16.2248	190.8418	24.8675	194.7933	101.3310	0.01
2007-03-27-1	13	1	16	16	187.5700	4.4339	179.5899	6.3662	186.8455	6.7636	204.2783	82.6860	0.01
2007-03-27-1	19	1	16	16	183.9616	4.7369	177.5880	9.0639	187.5362	12.5737	199.0751	75.2832	0.01
2007-03-27-1	12	14	8	8	96.1896	4.3270	91.5804	4.8478	101.0533	17.7677	95.2217	17.7215	0.05
2007-03-27-1	15	3	8	8	95.1950	5.3065	91.3280	4.4837	93.2359	14.0831	88.5928	12.4025	0.05
2007-03-27-1	18	1	8	16	94.8361	3.7355	89.5495	3.0290	90.1983	22.1848	88.7863	9.0227	0.05
2007-03-27-1	16	4	16	8	197.3781	39.0844	193.1085	45.8196	195.5270	46.9557	186.4053	44.2285	0.05
2007-03-27-1	13	1	16	16	194.0506	25.0185	189.7589	23.9145	178.5596	67.2600	190.7872	65.0006	0.05
2007-03-27-1	19	1	16	16	186.5073	32.0119	190.8256	38.0829	203.9420	64.1391	184.3914	43.0417	0.05
2007-03-27-1	12	14	8	8	94.5005	2.8645	93.2330	2.1564	95.5424	2.6671	95.3772	13.4551	0.005
2007-03-27-1	15	3	8	8	92.8686	2.3528	93.2227	2.3750	95.1239	2.8145	91.0736	12.8431	0.005
2007-03-27-1	18	1	8	16	93.0456	2.4397	91.0165	1.9198	93.8690	3.7720	90.9366	9.8092	0.005
2007-03-27-1	16	4	16	8	187.7490	13.0390	181.1116	18.0992	187.6466	20.4428	193.0291	49.7787	0.005
2007-03-27-1	13	1	16	16	187.7835	4.8401	181.8834	8.5070	187.6460	11.0390	194.4735	46.1891	0.005
2007-03-27-1	19	1	16	16	184.7671	5.1021	179.2986	11.9251	187.4085	12.4266	188.0205	36.1368	0.005
2007-03-27-1	12	14	8	8	94.7936	2.1326	91.3919	2.9559	96.8125	3.1910	98.7712	13.8309	0.03
2007-03-27-1	15	3	8	8	93.6018	2.9703	91.3809	2.9376	95.2784	4.0617	88.7937	12.5310	0.03
2007-03-27-1	18	1	8	16	93.7072	2.2236	89.9076	1.8954	92.7694	2.5143	90.1791	6.9358	0.03
2007-03-27-1	16	4	16	8	186.3865	23.7111	175.2925	36.8532	204.9289	57.0840	865.6901	3362	0.03
2007-03-27-1	13	1	16	16	189.1402	10.2458	176.6318	12.0661	179.5674	29.3187	190.4885	77.1416	0.03
2007-03-27-1	19	1	16	16	181.7615	11.5814	176.2752	21.9566	194.4906	54.3591	663.5999	3005.6	0.03
2007-03-27-1	12	14	8	8	94.0879	4.7149	92.5146	3.4401	95.8544	4.8683	97.7143	13.7864	0.003
2007-03-27-1	15	3	8	8	93.0628	3.8827	92.6951	4.6349	95.0814	7.2103	91.1907	11.0257	0.003
2007-03-27-1	18	1	8	16	93.1510	3.6988	91.3022	2.2069	93.9356	7.4918	93.5573	12.6705	0.003
2007-03-27-1	16	4	16	8	187.8483	15.7831	187.3727	21.2780	192.2736	21.8737	192.5124	32.2670	0.003
2007-03-27-1	13	1	16	16	187.6004	6.4725	182.1707	10.7201	188.1135	15.3255	183.5461	32.7766	0.003
2007-03-27-1	19	1	16	16	184.7493	7.7341	179.3949	13.8056	189.4849	17.6730	185.2179	33.3242	0.003];

dataBBMRpooled = [2007-8-24-4	10	2	4	8	48.2340	0.3298
2007-8-24-4	11	1	4	16	48.1423	0.5167
2007-8-24-4	2	9	8	8	96.1639	0.7256
2007-8-24-4	6	3	8	8	96.2724	0.9022
2007-03-27-1	12	14	8	8	94.5460	0.9911
2007-03-27-1	15	3	8	8	93.8937	1.0489
2007-8-24-4	7	3	8	16	96.1807	1.0616
2007-03-27-1	18	1	8	16	92.9802	1.1292
2007-8-24-4	4	1	16	8	191.6366	5.5426
2007-03-27-1	16	4	16	8	185.3886	7.8191
2007-8-24-4	5	4	16	16	191.0379	3.2867
2007-8-24-4	3	23	16	16	190.9040	2.9670
2007-03-27-1	13	1	16	16	186.4877	4.150
2007-03-27-1	19	1	16	16	184.1372	4.6096];

dataBBMLpooled = [2007-03-27-1	  12	16	8	8	94.9404	1.1027
2007-03-27-1	  15	2	8	8	93.9375	1.1196
2007-03-27-1	  18	2	8	16	92.4751	1.3550
2007-03-27-1	  16 	3	16	8	186.4490	9.5937
2007-03-27-1	  13	10	16	16	184.4530	4.3094
2007-03-27-1	  19	2	16	16	183.7410	5.6574
2007-8-24-4	  10	3	4	8	48.1397	0.3876
2007-8-24-4	  11	2	4	16	48.0469	0.6419
2007-8-24-4	  2	19	8	8	96.3183	0.7243
2007-8-24-4	  6	1	8	8	96.3395	1.0061
2007-8-24-4	  7	4	8	16	95.9980	0.9759
2007-8-24-4	  4	2	16	8	193.1800	5.9245
2007-8-24-4	  5 	3	16	16	191.7102	4.5259
2007-8-24-4	  3	17	16	16	192.0214	3.7722];

dataDBMRpooled = [2007-03-27-1	  12	2	8	8	93.8828	1.12
2007-03-27-1	  15	4	8	8	92.8571	1.6063
2007-03-27-1	  18	4	8	16	91.7633	1.5360
2007-03-27-1	  16 	1	16	8	185.7355	13.7776
2007-03-27-1	  13	2	16	16	182.6887	7.4844
2007-03-27-1	  19	4	16	16	180.6527	7.7815
2007-8-24-4	  10	4	4	8	48.0679	0.3875
2007-8-24-4	  11	4	4	16	48.3558	0.4075
2007-8-24-4	  2	7	8	8	96.5656	0.9908
2007-8-24-4	  6	4	8	8	96.7472	1.1963
2007-8-24-4	  7	1	8	16	96.4822	1.0562
2007-8-24-4	  4	3	16	8	192.7140	8.9752
2007-8-24-4	  5 	2	16	16	193.3761	4.2224
2007-8-24-4	  3	22	16	16	195.6119	4.0181];

dataDBMLpooled = [2007-03-27-1	  12	5	8	8	94.5537	1.3928
2007-03-27-1	  15	1	8	8	93.9146	1.7587
2007-03-27-1	  18	3	8	16	91.7963	1.6824
2007-03-27-1	  16 	2	16	8	188.1497	11.5963
2007-03-27-1	  13	7	16	16	187.9491	6.3612
2007-03-27-1	  19	3	16	16	186.2898	7.3682
2007-8-24-4	  10	1	4	8	48.1533	0.4902
2007-8-24-4	  11	3	4	16	47.6916	0.4933
2007-8-24-4	  2	17	8	8	95.7417	1.1051
2007-8-24-4	  6	2	8	8	95.7620	0.9915
2007-8-24-4	  7	2	8	16	95.5261	0.9628
2007-8-24-4	  4	4	16	8	191.2540	6.4294
2007-8-24-4	  5 	1	16	16	188.3077	4.7301
2007-8-24-4	  3	2	16	16	189.0769	4.5898];


dataALLpooled = [2007-8-24-4	10	2	4	8	48.2340	0.3298
2007-8-24-4	11	1	4	16	48.1423	0.5167
2007-8-24-4	2	9	8	8	96.1639	0.7256
2007-8-24-4	6	3	8	8	96.2724	0.9022
2007-03-27-1	12	14	8	8	94.5460	0.9911
2007-03-27-1	15	3	8	8	93.8937	1.0489
2007-8-24-4	7	3	8	16	96.1807	1.0616
2007-03-27-1	18	1	8	16	92.9802	1.1292
2007-8-24-4	4	1	16	8	191.6366	5.5426
2007-03-27-1	16	4	16	8	185.3886	7.8191
2007-8-24-4	5	4	16	16	191.0379	3.2867
2007-8-24-4	3	23	16	16	190.9040	2.9670
2007-03-27-1	13	1	16	16	186.4877	4.150
2007-03-27-1	19	1	16	16	184.1372	4.6096
2007-03-27-1	  12	2	8	8	93.8828	1.12
2007-03-27-1	  15	4	8	8	92.8571	1.6063
2007-03-27-1	  18	4	8	16	91.7633	1.5360
2007-03-27-1	  16 	1	16	8	185.7355	13.7776
2007-03-27-1	  13	2	16	16	182.6887	7.4844
2007-03-27-1	  19	4	16	16	180.6527	7.7815
2007-8-24-4	  10	4	4	8	48.0679	0.3875
2007-8-24-4	  11	4	4	16	48.3558	0.4075
2007-8-24-4	  2	7	8	8	96.5656	0.9908
2007-8-24-4	  6	4	8	8	96.7472	1.1963
2007-8-24-4	  7	1	8	16	96.4822	1.0562
2007-8-24-4	  4	3	16	8	192.7140	8.9752
2007-8-24-4	  5 	2	16	16	193.3761	4.2224
2007-8-24-4	  3	22	16	16	195.6119	4.0181
2007-03-27-1	  12	16	8	8	94.9404	1.1027
2007-03-27-1	  15	2	8	8	93.9375	1.1196
2007-03-27-1	  18	2	8	16	92.4751	1.3550
2007-03-27-1	  16 	3	16	8	186.4490	9.5937
2007-03-27-1	  13	10	16	16	184.4530	4.3094
2007-03-27-1	  19	2	16	16	183.7410	5.6574
2007-8-24-4	  10	3	4	8	48.1397	0.3876
2007-8-24-4	  11	2	4	16	48.0469	0.6419
2007-8-24-4	  2	19	8	8	96.3183	0.7243
2007-8-24-4	  6	1	8	8	96.3395	1.0061
2007-8-24-4	  7	4	8	16	95.9980	0.9759
2007-8-24-4	  4	2	16	8	193.1800	5.9245
2007-8-24-4	  5 	3	16	16	191.7102	4.5259
2007-8-24-4	  3	17	16	16	192.0214	3.7722
2007-03-27-1	  12	5	8	8	94.5537	1.3928
2007-03-27-1	  15	1	8	8	93.9146	1.7587
2007-03-27-1	  18	3	8	16	91.7963	1.6824
2007-03-27-1	  16 	2	16	8	188.1497	11.5963
2007-03-27-1	  13	7	16	16	187.9491	6.3612
2007-03-27-1	  19	3	16	16	186.2898	7.3682
2007-8-24-4	  10	1	4	8	48.1533	0.4902
2007-8-24-4	  11	3	4	16	47.6916	0.4933
2007-8-24-4	  2	17	8	8	95.7417	1.1051
2007-8-24-4	  6	2	8	8	95.7620	0.9915
2007-8-24-4	  7	2	8	16	95.5261	0.9628
2007-8-24-4	  4	4	16	8	191.2540	6.4294
2007-8-24-4	  5 	1	16	16	188.3077	4.7301
2007-8-24-4	  3	2	16	16	189.0769	4.5898];

if strcmp('BBMR', request)
    data = dataBBMR;
elseif strcmp('DBMR', request)
    data = dataDBMR;
elseif strcmp('BBML', request)
    data = dataBBML;
elseif strcmp('DBML', request)
    data = dataDBML;
elseif strcmp('ALL', request);
    data = [dataBBMR;dataDBMR; dataBBML; dataDBML];
elseif strcmp('tau', request);
    data = dataTau;
elseif strcmp('BBMRpooled', request);
    data = dataBBMRpooled;
elseif strcmp('BBMLpooled', request);
    data = dataBBMLpooled;
elseif strcmp('DBMRpooled', request);
    data = dataDBMRpooled;
elseif strcmp('DBMLpooled', request);
    data = dataDBMLpooled;
elseif strcmp('ALLpooled', request);
    data = dataALLpooled;  
end

    