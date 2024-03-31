function freeGB = freeDiskSpaceGB(path)

free = java.io.File(path).getFreeSpace();
freeGB = free/(1024*1024*1024);

end




