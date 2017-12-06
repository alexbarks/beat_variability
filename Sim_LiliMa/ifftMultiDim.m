%Takes in multi dimensional kspace array and converts to image space (if
%first two dimmensions are 2D image
function ksMx = ifftMultiDim(imgMx)
%4D array 
    [sy, sx, nPhases, nCoils] = size(imgMx); 
    ksMx = zeros(size(imgMx)); 
    
    
    for phaseNum = 1:nPhases
        for coilNum = 1:nCoils          
                ksMx(:,:, phaseNum, coilNum) = ifft2(fftshift(imgMx(:,:,phaseNum, coilNum))); 
        end
    end

end
