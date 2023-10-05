function [freq,Photons,PhotonMrk]=RetrievePhotons(LB,UB,cBlink,bin,PhotonTime)
    BinIndex=find((cBlink>=LB)&(cBlink<=UB));%BinIndex: the indices of cBlink, at which cBlink = l
    freq=length(BinIndex);%frequency: how many bins have cBlink in [LB:UB]
    PhotonMrk=ismember(bin,BinIndex);
    Photons=PhotonTime(PhotonMrk);%Photons: PhotonTime values corresponding to PhotonIndex
end