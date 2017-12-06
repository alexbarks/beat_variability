%Takes in multi dimensional kspace array and converts to image space (if
%first two dimmensions are 2D image
function imgMx = fftMultiDim(ksMx)
%4D array 
    [sy, sx, nPhases, nSlices, nCoils] = size(ksMx); 
    imgMx = zeros(size(ksMx)); 
    
    
    for phaseNum = 1:nPhases
        for sliceNum = 1:nSlices
            for coilNum = 1:nCoils
                imgMx(:,:, phaseNum, sliceNum, coilNum) = fftshift(fft2(ksMx(:,:,phaseNum, sliceNum, coilNum)));
            end
        end
    end

end
