function stacks = clear_tform(stacks, baseindex, goalindex)

base = get_stack(stacks, baseindex);
goalindex = parse_stack_index(goalindex);
base.tforms{goalindex{:}}     = [];
base.tforms_inv{goalindex{:}} = [];

stacks = set_stack(stacks, base, baseindex);

% Must rebuild routes to reflect deleted tform
stacks = build_tform_routes(stacks);