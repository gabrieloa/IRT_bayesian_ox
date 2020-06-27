

SBPtheta(const Resp) // Score Bruto Padronizado
{
decl media, sd, SB, SBP;

SB=sumc(Resp);

media = meanr(SB);
sd = sqrt(varr(SB));


SBP =  (SB - media)/ sd;

return (SBP);

}


