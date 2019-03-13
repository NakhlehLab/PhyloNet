package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.sitemodel;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.substitution.Frequencies;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.substitution.GTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.substitution.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by wendingqiao on 5/4/16.
 * Defines (1) mutation rate
 * (2) gamma distributed rates across sites (optional)
 * (3) proportion of the sites invariant (optional)
 * Currently assume no gamma/ no invariant
 */
public class SiteModel extends SiteModelInterface.Base {

    protected double _muParameter = 1.0;
//    private double _muLowerBound = 0; // not be used if fixed
//    private double _muUpperBound = Double.MAX_VALUE; // not be used if fixed
    protected int _gammaCategoryCount = 0;
    protected double _shapeParameter = Double.NaN; // shape parameter of gamma distribution. Ignored if gammaCategoryCount < 2
//    private double _shapeLowerBound = 1.0E-3;
//    private double _shapeUpperBound = 1.0E3;
    protected double _proportionInvariant = 0; // proportion of sites that is invariant, within [0, 1]
//    private double _proportionLowerBound = 0; // not be used if _proportionInvariant = 0
//    private double _proportionUpperBound = 1; // not be used if _proportionInvariant = 0

    protected boolean _ratesKnown;
    protected int _categoryCount = 0;
    protected double[] _categoryRates;
    protected double[] _categoryProportions;

    public SiteModel(SubstitutionModel model) {
        super(model);
        refresh();
    }

    public SiteModel(SubstitutionModel model, double mu, int gamma, double shape, double proportion) {
        super(model);
        _muParameter = mu;
        _gammaCategoryCount = gamma;
        _shapeParameter = shape;
        _proportionInvariant = proportion;
        if (_proportionInvariant < 0 || _proportionInvariant > 1) {
            throw new IllegalArgumentException("proportion invariant should be between 0 and 1");
        }
        refresh();
    }

    @Override
    protected void refresh() {
        if (!Double.isNaN(_shapeParameter)) {
            _categoryCount = _gammaCategoryCount;
            if (_categoryCount < 1) {
                if (_categoryCount < 0) {
                    System.err.println("SiteModel: Invalid category count (" + _categoryCount + ")! Set to 1");
                }
                _categoryCount = 1;
            }
        } else {
            _categoryCount = 1;
        }

        if (_proportionInvariant > 0) {
            if (_proportionInvariant >= 1.0) {
                throw new RuntimeException("Proportion invariant should be in bewteen 0 and 1");
            }
            if (_hasPropInvariantCategory) {
                _categoryCount += 1;
            }
        }
        _categoryRates = new double[_categoryCount];
        _categoryProportions = new double[_categoryCount];
        calculateCategoryRates(null, null);
    }


    @Override
    public double getProportionInvariant() {
        return _proportionInvariant;
    }

    @Override
    public boolean integrateAcrossCategories() {
        return true;
    }

    @Override
    public int getCategoryCount() {
        return _categoryCount;
    }

    @Override
    public int getCategoryOfSite(final int site, UltrametricTree tree, final TNode node) {
        throw new IllegalArgumentException("Integrating across categories");
    }

    @Override
    public double getRateForCategory(final int category, UltrametricTree tree, final TNode node) {
        synchronized (this) {
            if (!_ratesKnown) {
                calculateCategoryRates(tree, node);
            }
        }
        return _categoryRates[category] * _muParameter;
    }


    /**
     * return category rates
     *
     * @param node rates to which the rates apply. Typically, the rates will be uniform
     *             throughout the tree and the node argument is ignored.
     */
    @Override
    public double[] getCategoryRates(UltrametricTree tree, final TNode node) {
        synchronized (this) {
            if (!_ratesKnown) {
                calculateCategoryRates(tree, node);
            }
        }
        final double[] rates = new double[_categoryRates.length];
        for (int i = 0; i < rates.length; i++) {
            rates[i] = _categoryRates[i] * _muParameter;
        }

        return rates;
    }

    /**
     * Get the expected proportion of sites in this category.
     *
     * @param category the category number
     * @param node     node to which the proportions apply. Typically, proportions
     *                 will be uniform throughout the tree and this argument is ignored.
     * @return the proportion.
     */
    @Override
    public double getProportionForCategory(final int category, UltrametricTree tree, final TNode node) {
        synchronized (this) {
            if (!_ratesKnown) {
                calculateCategoryRates(tree, node);
            }
        }
        return _categoryProportions[category];
    }

    /**
     * Get an array of the expected proportion of sites in this category.
     *
     * @return an array of the proportion.
     */
    @Override
    public double[] getCategoryProportions(UltrametricTree tree, final TNode node) {
        synchronized (this) {
            if (!_ratesKnown) {
                calculateCategoryRates(tree, node);
            }
        }
        return _categoryProportions;
    }

    /**
     * discretisation of gamma distribution with equal proportions in each
     * category
     * @param node
     */
    protected void calculateCategoryRates(UltrametricTree tree, final TNode node) {
        double propVariable = 1.0;
        int cat = 0;

        if (_proportionInvariant > 0) {
            if (_hasPropInvariantCategory) {
                _categoryRates[0] = 0.0;
                _categoryProportions[0] = _proportionInvariant;
            }
            propVariable = 1.0 - _proportionInvariant;
            if (_hasPropInvariantCategory) {
                cat = 1;
            }
        }

        if (!Double.isNaN(_shapeParameter)) {
            // TODO
        } else {
            _categoryRates[cat] = 1.0 / propVariable;
            _categoryProportions[cat] = propVariable;
        }
        _ratesKnown = true;
    }


    @Override
    public double propose() {
        return 0;
    }

    @Override
    public void undo() {}

    @Override
    public void accept() {}

    @Override
    public void reject() {}

    @Override
    public double logDensity() {
        return 0;
    }

    @Override
    public boolean mayViolate() {
        return false;
    }

    @Override
    public boolean isValid() {
        return true;
    }

    // test
    public static void main(String[] args) {
        Map<String, String> omap = new HashMap<>();
        omap.put("L", "GTCATCCGACTCGTCTGACCCCGGCGCGGATCTGGAACCCCTGCCGTGACTCGGCGCGCAGGCGGGCGCATCGCCAAACCTAGCCCAGTGCAGGTCGTGGATGTAACGCTAGCGCTATTTGTGCATCAAAAGATAGCCGCTCGACGTCGCGCCCTTGTGGCCTTACCATCAGTGCGGCCCATCGGCCAGACCTCTGGAAGACTGAATGACAATATTTCGAAGGCGGAAGCCCGAAATGGCGTATGTGTCTGGAACAGAGGAGGGGTAGACTGCCCGTATGCACAGAACGCGACTGTGCGCTAGGCAAGGCTTGACCTAGCCGAATCGGTGCAGACTTTATTGCGACCGCAACGCTTAAGCCCCAGTGCGTCATGCCCCGTAAGTGTATGGGGCACTGCCCTTCGGCGCCCATCGTTGCCGCCTGGACGGCGACTCTCCTCATGGATTCCGCGCAACGGCCGGCTACCGCGAACGAGCCGTCATACACGTTTCGTACAGTA");
        omap.put("R", "GTTCACCGACTCGTCTGATCCCGACGCGGATCTGGGACCCATGCCATGGCTCAAGGCGCAGGCGGACGCATCGCCAAACCTAGCCCAGTGTAGGTCGTGGATGTCACACTAGGGCTATTTGTGCATCAGAAGATAGCTGCTCGACCTAGCGCCCTTGTTGCCCTGCCATCGGTGCGACCCATCGACCAGATCTCTGGAAGGCTGAACGATCATCTTGCGAAGGCGGAAGCCCGAAGTAACGTAGGTATCTGTAGTAGAAAAAGGGTAGACTACCCGTGTGCATAGAGTGCGACTGCGCGCTAGGCAAAGCTTGAACTAGCCGAATCGGTGCAGACCTTGTTGCGACCGCAATGCTTATGCCCCTGAGCGTCATGCGTCGTAAGTGTGTAAGGCACTGCGCTTCGGCGCCCCTCGCTGCCGCCTGGACGGCGCCTCTCCTCATGAATTCCGCGCAACGGCCGACTACCGCGAACGAGCCGCCATACACTTTTCGCGCATTG");
        omap.put("G", "GTCCTCCGACTCGTCTGATCCCGACGCGAACCCGGCGCGCCGGCCATGACTCGGGGCGCAGGCGGGCGCACCGTCAAGCCTAGCCCAGGGCAGATCGTGGATATCACGCTAGCGCTATTTGTGTATCAAAGCATAACTGCTCGGCGTAGCACCCTCGTGGCCTCTTCATCAGCGCGGCCCGTCGACCAGATTTCTGGAAGGCTGGATGACCATTTTGCGAGGGCGGGAGTCCGAAGTAGCGTAAGTGTCCGGAGCAGAGAAGGAGTAAACTGCCGGTGTGCTTAGAGCGCGACTGTGCGCCAGGCAAAACTTGAAATAGCCGAATCGATGCACACTCCGTTGCGACCGCAGTGCTTATGCCCCTGTGCGTCCTGCATCGTAAGGGTATAGGGCGCGGCGCTTCAGCGCCCGTCGTTGCCGCCTGGACGACGCCCCTTCTCACGGATGCCGCCCAACGGCCCGCTACCCCGAACGGTCAGCCAGACACTTTTCGTACATTA");
        omap.put("C", "GCCCTCCGACTCGTCTGATCCCGACGCGGACCCGACGCACCGGCCATGACTCGGGGCGCAGGCGGGCGCACCGTCAAGCCTAGCCCAGGGCAGACCGTGGATATCGCGCTAGCGCTATTTGTGCATCAAAACGTAACTGCTCGGCGTAGCACCCTTGTGGCCGCTCCGTCAGCGCGGCCCGTCGACCAGATTTCTGGAAGGCTGGATGACCATCTTGCGGGGGCGGGAGTCCGCAGTAGCGTAAGTGTCCGGAGCAGAGAAGGAGTAAACTGCCGGTATTCATAAAGTGCAACTGTGCGCCAGGCAAAACTTGAAATAGCCGATTCGATGCACACCTCGTTGCAACCGAAGTGCCTATGCCCCTGTGCGTCCTGCATCGTAAGGGTATAGGGCGCGGCGCTTCAGCGCCCACCGTTGCCGTCTGAACGGCGCCTCTTCTCACGGATGCCGCCCAACGGCCGGCTACCCCGAACGGTCAGCCAGGCACTTTTCGTACATTA");
        omap.put("Q", "GTCCTTCGACTCGTCTGATCCCGACGCGGACCCGGCGCCCCTGCCATGACTCGGGGCGCAGGCAGGCGCACCGTCAAATCGAGCCCAGTGCAGGTCGTGGATGTCACGCTAGCGCTATCTGTGCATCAAGGGGTAACTGTTCGCCGTAGTGCCCTTGTGGCCTCGCCATCAGTGCGGCCCATCGACCAGGTCTCTGGAAGGCTGGATAACCACTTTGCGAAGACGGGAGTCCGAAGTGGCGTATGTGTCCAGAGCAGAGAAGGGATAAACTGCCCGAATGCGTAGAGCGCGACCGTGCGCTAGGCAAAGTCTGAGCTAGCCGAATCGGCGCGCACTTTATTGCGGCCGCAGTGCCTTGGTCCCTGTGCGTCCCGCAACGTAAGGGTATGGGGCACCGCGCTTCGGCGCCTATCCTTGCCGCCTGGACGGCGCTTCTTCTCATGGATTCCACGCAACGACCGGCACCTGCGAACGGTCCGCCAGACACTTTTCGTACACTA");
        omap.put("A", "GTCCTCCGACTCATCTGATCCCGACGCGGACCCGGCGCCCCGGCCATGACTCGGGGCGCAGGCGGGCGCACCGTCAAGCCTTGCCCGGTGCAGGTCGTGGATGTCACGCTAGCGCTATCTGCGCATCAAAAGGTAACTGTTCGACGTGGTGCCCTTGTGGCCTCGCCATCAGTGCGGCCCATCGACCAGGTCTCTGGAAGGCTGTATAACCACTTTGCGACGACGGGAGTCCGAAGTGGCGTAAGTGTCCGGAGCAGAGAAGGGATAAACTGCCTGAATGCGCAGAGCGCGACCGTGCGCTAGGCAAAGTTTGAGCTAGCCGAATCGGCGCGCACTTTATTGCGGCCGCAGTGTGTTGGTCCCTGTGCGTCCCGCAACGTAAGGGTATGGGGCACCGCGCTTCGGCGCCTATCGTTGCCGCCTGTACGGCGATTCTTCTCATGAATCCCACACATCGACCGGCACCTGCGAACGGTCCGCCAGACACTTTTCGTACACTA");

        MarkerSeq aln = new MarkerSeq(omap);

        double[] base = {0.2112, 0.2888, 0.2896, 0.2104};
        Frequencies freq = new Frequencies(aln, base);

        double[] trans = {0.2173/0.2070, 0.9798/0.2070, 0.2575/0.2070, 0.1038/0.2070, 1/0.2070, 1.0};
        SubstitutionModel subst = new GTR(freq, trans);

        SiteModel site = new SiteModel(subst);
        System.out.println(site.getCategoryCount() == 1);
        System.out.println(Arrays.toString(site.getCategoryProportions(null, null)).equals("[1.0]"));
    }
}