

SBPb(const Resp) // Score Bruto Padronizado
{
decl media, sd, SB, SBP;

SB=sumr(Resp);

media = meanc(SB);
sd = sqrt(varc(SB));


SBP =  (SB - media)/ sd;

return (-SBP);

}


