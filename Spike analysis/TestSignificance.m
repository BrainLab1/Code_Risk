

function [p , sampleSide] = TestSignificance(allSamples, oneSample, alpha, testType)


x1 = sum(allSamples <= oneSample)/length(allSamples);
x2 = sum(allSamples > oneSample)/length(allSamples);

switch testType
    
    case 'two side'
        if x1 <= alpha/2
            sampleSide = 'Left';
            p = x1;
            elseif x2 <= alpha/2
                sampleSide = 'Right';
                p = x2;
        else
                sampleSide = 'Center';
                p = NaN;
        end
                    
end % end of testType


return


