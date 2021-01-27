function [RegulatedModel, F_KO, TFs] = PRIME (model, Regulator, Target, Magnitude, GeneExpressionIDs, GeneExpression, overallExpressionid, minOverallExpression, maxOverallExpression)
%%
%PRIME is an algorithm used to integrate regulatory network onto metabolic networks
%It predicts the phenotype of TF and gene knockouts
%It also gives the perturbed models as output
%PRIME calculates the magnitude of TF influences on the genes and their associated reactions
%
%USAGE:
%   [RegulatedModel, F_KO, TFs] = PRIME (model, Regulator, Target,Magnitude, GeneExpressionIDs, GeneExpression, overallExpressionid, minOverallExpression, maxOverallExpression);
%
%INPUT:
%   model:                  Metabolic network as COBRA structure (with expression data applied using iMAT/GIMME) e.g., iEK1011
%   Regulator and Target:   Cell arrays listing the TF-target gene interactions. Regulators = {list of TFs}; Targets = {list of gene targets};
%                           Regulators and Targets have the same number of rows; each row describes a separate interaction.
%   Magnitude:              Weights (Beta) from INFERELATOR run of EGRIN. Magnitude has the same number of rowss as Regulators and Targets
%   GeneExpressionIDs       Gene ID of GE data (present study) e.g., isoniazid treatment
%   Gene Expression         GE data corresponding to the GeneExpressionIDs (present study)
%   overallExpressionid     Gene ID of GE data (Gene expression from various conditions)
%   minOverallExpression    Min value of GE data from various conditions
%   maxOverallExpression    Max value of GE data from various conditions
%
%OUTPUT:
%   RegulatedModel  TF pertubed model (Models of TF KOs)
%   F_KO:           Growth rate of Knockouts (Objective function of the model)
%   TFs:            List of TFs being knocked out or over expressed
%
% .. Authors:
%   - Selva Rupa Christinal Immanuel, ISB, Seattle, USA. 05/31/2019 Original code
%% COBRA toolbox intialization
xx=waitbar(0,'PRIME in progress ...');
%disp('COBRA toolbox initialization...')
%initCobraToolbox
%changeCobraSolver('glpk');
%% INPUT processing for PRIME run
disp('processing inputs...')
threshold = 10^(-6);
[uniqueTFs] = unique(Regulator);
model = buildRxnGeneMat(model);
waitbar(0.15,xx,'PRIME in progress ...');
disp('running Flux Variability Analysis...')
[minFlux, maxFlux] = fluxVariability(model); %gives the best flux span for the condition
%minFlux = model.lb; %default minFlux
%maxFlux = model.ub; %default maxFlux
minFlux(abs(minFlux) < threshold) = 0;
maxFlux(abs(maxFlux) < threshold) = 0;
Magnitude = (Magnitude-min(Magnitude))/(max(Magnitude)-min(Magnitude)); %rank the magnitude value to compare the models accrosss conditions
%%
disp('initializing PRIME...')
waitbar(0.33,xx,'PRIME in progress ...');
%objective = [find(model.c) obj_frac]; 
RxnNumber = size(model.S,2);
%parsedGPR = GPRparser(model);
%% Running PRIME
%Apply Refine value of each TF perturbation to its corresponding reactions 
TFs = uniqueTFs;
genelist = ismember(Target,model.genes);
nGenes = length(model.genes);
nRxns = length(model.rxns);
disp('calculating TF influences...')
fba = optimizeCbModel(model);
WTGr = fba.f;
fbaFlux = fba.x;
for a = 1:RxnNumber
    Refine = ones(a,1);
end
waitbar(0.5,xx,'PRIME in progress ...'); 
%litevidence = logical(Litevidence);
%Probability_gene(litevidence) = Probability(litevidence);
[ReactionPosition,genePosition] = find(model.rxnGeneMat);
sumgeneBetas = ones(nGenes,1);

for n = 1:length(model.genes)
    genePos = ismember(Target,model.genes(n));
    geneMagnitude = Magnitude(genePos);
    %geneMagnitude = (geneMagnitude - min(geneMagnitude))/(max(geneMagnitude)-min(geneMagnitude));
    sumgeneBetas(n,1) = sum(abs(geneMagnitude));
    %sumgeneBetas(n,1) = sum(geneMagnitude);
end
waitbar(0.75,xx,'PRIME in progress ...');    
disp('mapping ReFInE to maxFlux of Reactions...')
disp('running Flux Balance Ananlysis on the perturbed model')
for m = 1:length(TFs)
    TFposition = ismember(Regulator,TFs(m));
    Targets = Target(TFposition); 
    %TargetProbability = Probability_gene(TFposition);
    TargetMagnitude = Magnitude(TFposition);
    %TargetMagnitude = (TargetMagnitude - min(TargetMagnitude))/(max(TargetMagnitude)-min(TargetMagnitude));
    nReactions = length(ReactionPosition);
    %Likelihood = ones(nReactions,1);
    calculatedBetas = ones(nReactions,1);
    Factor = ones(nReactions,1);
    FactorMultiplied = ones(nReactions,1);
    GEdata = mean(GeneExpression,2);
    TFpositionGE = ismember(GeneExpressionIDs, TFs(m));
    TFinGE = find(TFpositionGE,1);
            if ~isempty(TFinGE)
                     TFactivity(m,1)= GEdata(TFpositionGE);
                else TFactivity(m,1)= min(GEdata);
            end
                TFpositionOverallGE = ismember(overallExpressionid, TFs(m));
                TFinOverall = find(TFpositionOverallGE,1);
                 if ~isempty(TFinOverall)
                    TFminOverall = minOverallExpression(TFpositionOverallGE);
                    TFmaxOverall = maxOverallExpression(TFpositionOverallGE);
                 else TFminOverall && TFmaxOverall == 0;
                 end
                        if TFactivity (m,1) == 0
                                 TFactivityCondition(m,1) = TFactivity(m,1);
                            else TFactivityCondition(m,1) = (TFactivity(m,1)-TFminOverall)/(TFmaxOverall-TFminOverall);
                        end
            for i =1:length(ReactionPosition)
                gene = model.genes(genePosition(i));
                geneID = find(ismember(Targets,gene));
                geneIDmodel = find(ismember(model.genes,gene));
                if ~isempty (geneID)
                    %Likelihood(i) = TargetProbability(geneID);
                    Beta = TargetMagnitude(geneID);
                    calculatedBetas(i)= 1-((Beta/sumgeneBetas(geneIDmodel))*TFactivityCondition(m,1)); %ReFInE formula
                    if calculatedBetas(i) == 0
                       calculatedBetas(i) = 1-(Beta*TFactivityCondition(m,1));
                    end
                end
            end
                    Factor = calculatedBetas;
                    FactorMultiplied = FactorMultiplied.*Factor;%all Beta associated to the reaction
                    GPRreactions = unique(ReactionPosition);
                    Reactions = model.rxns(GPRreactions);
                    for  k = 1:length(GPRreactions)
                        GPRreactionPosition = ismember(ReactionPosition,(GPRreactions(k)));
                        ReactionIndex = find(GPRreactionPosition);
                        FinalFactor = min(FactorMultiplied(ReactionIndex)); %consider the min of all genes involved in the reaction
                        ModelReactionIndex = ismember(model.rxns,Reactions(k));
                        ModelReactionPosition = find(ModelReactionIndex);
                        Refine(ModelReactionPosition) = FinalFactor; %vector of Refine; assigned to the reaction position
                    end
    RegulatedModel{m,1} = model;
    inhibit = RegulatedModel{m,1}.c;
    RegulatedModel{m,1}.c = inhibit.*min(abs(Refine)); % new weight factor
    RegulatedModel{m,1}.ub = maxFlux.*Refine;%max flux perturbation using Refine
    RegulatedModel{m,1}.lb = minFlux.*Refine;%min flux perturbation using Refine
    ubound = RegulatedModel{m,1}.ub;
    lbound = RegulatedModel{m,1}.lb;
                    for x = 1:RxnNumber
                        ubound(x,1) = max(ubound(x,1),fbaFlux(x,1));%maintain stoichiometry consistency to be lb<ub
                        lbound(x,1) = min(lbound(x,1),fbaFlux(x,1));%maintain stoichiometry consistency to be lb<ub
                    end
    RegulatedModel{m,1}.ub = ubound;
    RegulatedModel{m,1}.lb = lbound;
    FBAsoln = optimizeCbModel(RegulatedModel{m});
    F_KO(m,1) = FBAsoln.f;
    RefineValue(:,m)=Refine;
end
waitbar(0.9,xx,'PRIME in progress ...');
close(xx)
end                     
     

