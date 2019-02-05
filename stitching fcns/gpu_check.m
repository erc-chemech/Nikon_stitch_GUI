% Determine GPU capabilities. If it is sufficient, use GPU device.
function [dev,flag4]=gpu_check()
    if license('test','Distrib_Computing_Toolbox')==0
        flag4=0;
        dev=nan;
        return
    end
    if ~(parallel.gpu.GPUDevice.isAvailable)
        fprintf(['\n\t**GPU does not have a compute capability of 1.3 or '...
                 'greater.']);
        dev=[];
        flag4=false;
    else
        dev = gpuDevice;
        fprintf(...
        'GPU detected (%s, %d multiprocessors, Compute Capability %s)',...
        dev.Name, dev.MultiprocessorCount, dev.ComputeCapability);
        flag4=true;
    end
end