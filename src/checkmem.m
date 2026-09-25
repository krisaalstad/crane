workspace_info = whos;
total_memory_bytes = sum([workspace_info.bytes]);
MB = total_memory_bytes / (1024^2);
fprintf('Total workspace memory: %.2f MB\n', MB);


