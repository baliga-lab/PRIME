function [FinalList, ListsRxn, ListsGn, BetaFinal] = PRIMEanalysis (model, regulator, target, TF, Magnitude)
%INPUT:
%   model:                  Metabolic network as COBRA structure (with expression data applied using iMAT/GIMME)
%   regulator:              Cell array listing regulators from PRIME
%   target:                 Cell array listing gene targets of regulators
%   TF:                     TF or list of TFs for which the analysis is being done
%   Magnitude               Weights from EGRIN
%OUTPUT:
%   Final List: Table of TF, gene target, reactions and pathways
%   ListRxn:  List of reactions associated with the TF or TFs
%   ListGn: List of gene targets for the TF in the metabolic model
%   BetaFinal: Table of TF, Gene and the weights associated with the interaction

model = buildRxnGeneMat(model);
GeneList = model.genes;
PathwayList = model.subSystems;
ReactionList = model.rxns;
ReactionNameList = model.rxnNames;
[ReactionPosition,genePosition] = find(model.rxnGeneMat);
ListsRxn = {};
ListsGn = {};
FinalList = {};
BetaFinal = [];
TFname = {};
TF = unique(TF);
for m = 1:length(TF)
TFposition = ismember(regulator,TF(m));
Targets = target(TFposition);
GeneIDs = ismember(GeneList,Targets);
fGeneIDs = find(GeneIDs);
Genes = GeneList(GeneIDs);
Geneposition = ismember(Targets, Genes);
Gene = Targets(Geneposition);
regulatorM = regulator(TFposition);
regulatorList = regulatorM(Geneposition);
BetaM = Magnitude(TFposition);
betaposition = ismember(Targets,Genes);
Beta = BetaM(betaposition);
ReactionID = ReactionPosition(ismember(genePosition,fGeneIDs)); 
Reactions = ReactionList(ReactionID);
GeneID = genePosition(ismember(genePosition,fGeneIDs)); 
GenesFinal = GeneList(GeneID);
ReactionNames = ReactionNameList(ReactionID);
Pathways = PathwayList(ReactionID);
TFname = {};
TFnameGene = {};

    if ~isempty(Reactions)
      for i = 1: length(Reactions)
          for n = 1: length(GenesFinal)
          TFname(i,:) = TF(m);
          TFnameGene (n,:) = TF(m);
          end
      end
    TFRxnList = [TFname,Reactions,ReactionNames,Pathways];
    TFgeneList = [TFnameGene,GenesFinal];
    ListsRxn = [ListsRxn;TFRxnList];
    ListsGn= [ListsGn;TFgeneList];
    BetaF = [regulatorList,Gene,Beta];
    BetaFinal = [BetaFinal;BetaF];
    List = [TFname,GenesFinal,Reactions,ReactionNames,Pathways];
    FinalList = [FinalList;List];
    end
end
end