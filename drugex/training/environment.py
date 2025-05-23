import pandas as pd
from rdkit import Chem

from drugex.training.interfaces import Environment
from drugex.training.scorers.smiles import SmilesChecker


class DrugExEnvironment(Environment):
    """
    Original implementation of the environment scoring strategy for DrugEx v3.

    """

    def getScores(self, smiles, frags=None, no_multifrag_smiles=True):

        """
        This method is used to get the scores from the scorers and to check molecule validity and desireability.
        
        Parameters
        ----------
        smiles : list of str
            List of SMILES strings to score.
        frags : list of str, optional
            List of SMILES strings of fragments to check for.
        no_multifrag_smiles : bool, optional
            Whether to check for SMILES strings that contain more than one fragment.
        
        Returns
        -------
        scores : pd.DataFrame
            Dataframe with the scores from the scorers and the validity and desireability of the molecules.
        """
        mols = [Chem.MolFromSmiles(s) if s != '' else None for s in smiles]
        
        # Check molecule validity and accuracy
        scores = SmilesChecker.checkSmiles(smiles, frags=frags, no_multifrag_smiles=no_multifrag_smiles)

        # Get scores per objective from the scorers
        for scorer in self.scorers:
            scores.loc[:, scorer.getKey()] = scorer(mols)

        # Set all scores to 0.0 for invalid molecules
        scores.loc[scores['Valid'] == 0, self.getScorerKeys()] = 0.0

        # Check if the molecule is desirable
        undesire = (scores[self.getScorerKeys()] < self.thresholds)
        scores['Desired'] = (undesire.sum(axis=1) == 0).astype(int)

        # TODO: Maybe smiles and frags should be added to the dataframe as well?
        
        return scores

    def getUnmodifiedScores(self, smiles):

        """
        This method is used to get the scores from the scorers without any modifications.
        
        Parameters
        ----------
        smiles : list of str
            List of SMILES strings to score.
        
        Returns
        -------
        preds : pd.DataFrame
            Dataframe with the scores from the scorers.
        """
        
        preds = {}
        mols = [Chem.MolFromSmiles(s) for s in smiles]
        for scorer in self.scorers:
            score = scorer.getScores(mols)
            if isinstance(scorer.getKey(), list):
                for i, key in enumerate(scorer.getKey()):
                    preds[key] = score[:, i]
            else:
                preds[scorer.getKey()] = score
        preds = pd.DataFrame(preds)

        return preds