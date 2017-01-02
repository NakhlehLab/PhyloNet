package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.sitemodel;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.datatype.DataType;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.substitution.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by wendingqiao on 5/4/16.
 */
public interface SiteModelInterface {

    void setDataType(DataType dataType);

    // Base implementation of a site model with substitution model and rate categories
    abstract class Base extends StateNode implements SiteModelInterface {

        protected DataType m_dataType;
        protected double _proportionInvariant = 0.0;
        protected SubstitutionModel _substModel = null;

        public boolean _hasPropInvariantCategory = true;

        protected List<StateNode> _conditions = null;

        public Base(SubstitutionModel model) {
            _substModel = model;
        }

        /**
         * Specifies whether SiteModel should integrate over the different categories at
         * each site. If true, the SiteModel will calculate the likelihood of each site
         * for each category. If false it will assume that there is each site can have a
         * different category.
         */
        abstract public boolean integrateAcrossCategories();

        /**
         * Gets the number of categories of substitution processes
         */
        abstract public int getCategoryCount();

        /**
         * Gets the category of a particular site.
         */
        abstract public int getCategoryOfSite(int site, UltrametricTree tree, TNode node);

        /**
         * Gets the rate for a particular category.
         */
        abstract public double getRateForCategory(int category, UltrametricTree tree, TNode node);

        /**
         * Gets an array of the rates for all categories.
         */
        abstract public double[] getCategoryRates(UltrametricTree tree, TNode node);

        /**
         * Gets the expected proportion of sites in this category.
         */
        abstract public double getProportionForCategory(int category, UltrametricTree tree, TNode node);

        /**
         * Gets an array of the expected proportion of sites for all categories.
         */
        abstract public double[] getCategoryProportions(UltrametricTree tree, TNode node);

        public void setPropInvariantIsCategory(final boolean propInvariantIsCategory) {
            _hasPropInvariantCategory = propInvariantIsCategory;
            refresh();
        }

        /**
         * set up categories, reserve appropriately sized memory
         */
        protected void refresh() {}

        /**
         * Gets this site model's substitution model
         */
        public SubstitutionModel getSubstitutionModel() {
            return _substModel;
        }

        public List<StateNode> getConditions() {
            return _conditions;
        }

        public void addCondition(final StateNode stateNode) {
            if (_conditions == null) _conditions = new ArrayList<>();
            _conditions.add(stateNode);
        }

        @Override
        public void setDataType(final DataType dataType) {
            m_dataType = dataType;
        }

        public double getProportionInvariant() {
            return _proportionInvariant;
        }

    }

}
