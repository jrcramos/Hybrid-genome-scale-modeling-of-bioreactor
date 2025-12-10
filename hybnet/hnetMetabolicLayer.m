function layer=hnetMetabolicLayer(ninputs,numnodes,Si,Se,Isrev,inpl)%###############
if nargin<6
    inpl=-1;
end
layer.type = 'Metabolic';
layer.nx = ninputs;
layer.ny = numnodes;
layer.inpl=inpl;
layer.Si = Si; %intracellular stoichiometry
layer.Se = Se; %extracellular stoichiometry
layer.Ir = ~Isrev; %irreversible reactions
layer.nw = 0;
layer.nullsi = null(Si);    % v = nullsi * t
layer.senullsi = Se*layer.nullsi;   %re = senullsi * t
layer.senullsi_t = layer.senullsi';   
layer.initwfun = @metainitw;
layer.setwfun = @metasetw;
layer.getwfun = @metagetw;
layer.getgradfun = @metagetgrad;
layer.bpresetfun = @metaresetbp;
layer.updwfun=@metaupdweights;
layer.fffun = @metaff;
layer.bpfun = @metabp;
layer.initstatefun=@metainitstate;
assert(layer.nx==size(layer.nullsi,2),...
    'Number of inputs does not match the number of columns of Null(Si)');
assert(layer.ny==size(layer.Se,1),...
    'Number of outputs does not match the number of rows of Se');
end

function layer=metainitw(layer,dw)%########################################
end

function layer=metainitstate(layer)%########################################
end

function layer=metasetw(layer,w)%############################################
end

function w=metagetw(layer)%##################################################
w=[];
end

function [dlossdw,sum_dlossdw]=metagetgrad(layer)%###############################
dlossdw=[];
sum_dlossdw=[];
end

function layer=metaff(layer,x)%##############################################
layer.x=x;
layer.v = layer.nullsi*layer.x;
layer.y = layer.Se*layer.v;
layer.dlossdx=zeros(layer.nx,1); %exception  !!!!!!!!!!!!!!!!!!!
end

function layer=metabp(layer,dlossdy)%######################################
layer.dlossdw=[];
%layer.y=layer.senullsi*layer.x
%layer.dlossdx=layer.senullsi_t*dlossdy;
layer.dlossdx=layer.dlossdx+layer.senullsi_t*dlossdy; %exception
end

function layer=metaresetbp(layer)%#########################################
layer.dlossdx=zeros(layer.nx,1);
layer.dlossdy=zeros(layer.ny,1);
end

function layer=metaupdweights(layer,iter,alfa,beta1,beta2,eta,dropout)%#####
end