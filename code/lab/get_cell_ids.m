function cell_ids=get_cell_ids(dataset,cell_identifier);

    index=get_cell_indices(dataset,cell_identifier);

    cell_ids=dataset.cell_ids(index);