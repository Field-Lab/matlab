
% open the old model file to get some clusters out of it (this is an
% example, you could generate clusters in any way you want).
electrode = 57;
model = edu.ucsc.neurobiology.vision.io.ClusteringModelFile('/Volumes/Analysis/2012-08-09-3/data005-from-data002/data005-from-data002.model');
% get the clusters on electrode 1
clusters = model.getNeuronExtraction(electrode);
clusters.probability

%%
% model
% newModelObject.addElectrode(clusters.probability, clusters.means, clusters.covariances, electrode, nClusters)
% newModelObject.closeModel;

% plot the output verifying that everything worked
% projectionsPath = '/Volumes/Analysis/2012-08-09-3/data002/data002.prj';
% figure; plot_projections(projectionsPath, 'clusters', clusters, 'electrode', electrode);
