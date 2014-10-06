function [C1,C2] = getC1C2( FirstArtifact,SecondArtifact,FirstMovieAmplitude,SecondMovieAmplitude,Kropka)
% Ta funkcja wyznacza wartoœci wspó³czynnników C1 oraz C2

if Kropka ==1

    C1=FirstArtifact-((FirstArtifact-SecondArtifact)./(FirstMovieAmplitude-SecondMovieAmplitude)).*FirstMovieAmplitude;   %%  ??? czy .*
    C2=(FirstArtifact-SecondArtifact)./(FirstMovieAmplitude-SecondMovieAmplitude);

else
    
    C1=FirstArtifact-((FirstArtifact-SecondArtifact)/(FirstMovieAmplitude-SecondMovieAmplitude))*FirstMovieAmplitude;   %%  ??? czy .*
    C2=(FirstArtifact-SecondArtifact)/(FirstMovieAmplitude-SecondMovieAmplitude);

end
end

