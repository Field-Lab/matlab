function pth = get_stack_basepath(stack)

pth = '';
if isempty(stack) || ~isfield(stack, 'basepath') || isempty(stack.basepath)
    return;
end

if isa(stack.basepath, 'function_handle')
    pth = stack.basepath();
    return;
end

pth = stack.basepath;