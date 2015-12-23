function path = reconstructPath(camefrom, current)

if camefrom(current)
    path = reconstructPath(camefrom, camefrom(current));
    path = cat(2, current, path);
else
    path = current;
end

end