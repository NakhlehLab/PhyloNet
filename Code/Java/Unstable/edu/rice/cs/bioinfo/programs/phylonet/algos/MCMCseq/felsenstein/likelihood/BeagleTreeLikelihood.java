package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.likelihood;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.BeagleInfo;
import beagle.InstanceDetails;
import beagle.ResourceDetails;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.StateNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.branchRate.BranchRateModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.branchRate.StrictClockModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.sitemodel.SiteModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.EigenDecomposition;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.Frequencies;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.GTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.substitution.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by wendingqiao on 5/5/16.
 * Use Beagle library to calculate "Felsenstein likelihood".
 * Adapted from
 * https://github.com/CompEvol/beast2/blob/master/src/beast/evolution/likelihood/BeagleTreeLikelihood.java
 */
public class BeagleTreeLikelihood extends StateNode {

    private boolean _print = false;

    // --- distribution ---
    protected double _logP = Double.NaN;
    protected double _storedLogP = Double.NaN;

    // --- generic tree likelihood ---
    protected Alignment _alignment;
    protected UltrametricTree _geneTree;
    protected SiteModel _siteModel;
    protected BranchRateModel _branchRateModel;

    // --- tree likelihood ---
    protected SubstitutionModel _substitutionModel;

    // special cases
    protected boolean _useAmbiguities = false;
    protected boolean _useTipLikelihoods = false;
    protected boolean _useAscertainedSitePatterns = false;
    public static enum Scaling {none, always, _default};
    protected Scaling _scaling = Scaling._default;
    protected double _proportionInvariant = 0;
    protected List<Integer> _constantPattern = null;
    protected int _invariantCategory = -1;

    protected int _nStateCount;
    protected int _nNodeCount;
    protected int _categoryCount; // # rate categories

    private double [] _currentCategoryRates;
    private double [] _currentFreqs;
    private double [] _currentCategoryWeights;


    /**
     * flag to indicate the
     * - CLEAN=0: nothing needs to be recalculated for the node
     * - DIRTY=1: a node partial needs to be recalculated
     * - FILTHY=2: the indices for the node need to be recalculated
     */
    protected int _hasDirt;

    protected double[] _prevBranchLengths;
    protected double[] _storedBranchLengths;

    //store likelihoods for each pattern
    protected double[] _patternLogLikelihoods;
    // store partials
    protected double[] _fRootPartials;
    protected double[] _tipPartials;
    // store probabilities
    protected double[] _probabilities;
    protected int _matrixSize;

    // --- beagle tree likelihood ---
    protected Beagle _beagle; // the BEAGLE library instance
    // BEAGLE properties
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";
    // Which scheme to use if choice not specified
    public enum PartialsRescalingScheme {
        DEFAULT("default"), // whatever our current favourite default is
        NONE("none"),       // no scaling
        DYNAMIC("dynamic"), // rescale when needed and reuse scaling factors
        ALWAYS("always"),   // rescale every node, every site, every time - slow but safe
        DELAYED("delayed"), // postpone until first underflow then switch to 'always'
        AUTO("auto");       // BEAGLE automatic scaling - currently playing it safe with 'always'

        private final String text;
        PartialsRescalingScheme(String text) {
            this.text = text;
        }
        public String getText() {
            return text;
        }
        public static PartialsRescalingScheme parseFromString(String text) {
            for (PartialsRescalingScheme scheme : PartialsRescalingScheme.values()) {
                if (scheme.getText().compareToIgnoreCase(text) == 0)
                    return scheme;
            }
            return DEFAULT;
        }
    }
    private static final PartialsRescalingScheme DEFAULT_RESCALING_SCHEME = PartialsRescalingScheme.DYNAMIC;

    private static int _instanceCount = 0;
    private static List<Integer> _resourceOrder = null;
    private static List<Integer> _preferredOrder = null;
    private static List<Integer> _requiredOrder = null;
    private static List<String> _scalingOrder = null;

    private static final int RESCALE_FREQUENCY = 10000;
    private static final int RESCALE_TIMES = 1;


    // --- temporary variables ---
    private int _eigenCount;
    private int[][] _matrixUpdateIndices;
    private double[][] _branchLengths;
    private int[] _branchUpdateCount;
    private int[] _scaleBufferIndices;
    private int[] _storedScaleBufferIndices;

    private int[][] _operations;
    private int _operationListCount;
    private int[] _operationCount;

    // Flag to specify that the substitution model has changed
    protected boolean _updateSubstitutionModel;
    protected boolean _storedUpdateSubstitutionModel;
    // Flag to specify that the site model has changed
    protected boolean _updateSiteModel;
    protected boolean _storedUpdateSiteModel;

    protected int _tipCount;
    protected int _internalNodeCount;
    protected int _patternCount;

    private PartialsRescalingScheme _rescalingScheme = DEFAULT_RESCALING_SCHEME;
    private int _rescalingFrequency = RESCALE_FREQUENCY;
    protected boolean _useScaleFactors = false;
    private boolean _useAutoScaling = false;
    private boolean _recomputeScaleFactors = false;
    private boolean _everUnderflowed = false;
    private int _rescalingCount = 0;
    private int _rescalingCountInner = 0;

    // --- BEAGLE Buffers ---
    protected BufferIndexHelper _partialBufferHelper;
    private BufferIndexHelper _eigenBufferHelper;
    protected BufferIndexHelper _matrixBufferHelper;
    protected BufferIndexHelper _scaleBufferHelper;

    protected class BufferIndexHelper {

        private final int maxIndexValue;
        private final int minIndexValue;
        private final int offsetCount;
        private int[] indexOffsets;
        private int[] storedIndexOffsets;

        /**
         * @param maxIndexValue the number of possible input values for the index
         * @param minIndexValue the minimum index value to have the mirrored buffers
         */
        BufferIndexHelper(int maxIndexValue, int minIndexValue) {
            this.maxIndexValue = maxIndexValue;
            this.minIndexValue = minIndexValue;
            offsetCount = maxIndexValue - minIndexValue;
            indexOffsets = new int[offsetCount];
            storedIndexOffsets = new int[offsetCount];
        }

        public int getBufferCount() {
            return 2 * offsetCount + minIndexValue;
        }

        void flipOffset(int i) {
            if (i >= minIndexValue) {
                indexOffsets[i - minIndexValue] = offsetCount - indexOffsets[i - minIndexValue];
            }
        }

        int getOffsetIndex(int i) {
            if (i < minIndexValue) {
                return i;
            }
            return indexOffsets[i - minIndexValue] + i;
        }

        void getIndices(int[] outIndices) {
            for (int i = 0; i < maxIndexValue; i++) {
                outIndices[i] = getOffsetIndex(i);
            }
        }

        void storeState() {
            System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);

        }

        void restoreState() {
            int[] tmp = storedIndexOffsets;
            storedIndexOffsets = indexOffsets;
            indexOffsets = tmp;
        }
    }

    public BeagleTreeLikelihood(Alignment aln, UltrametricTree tree,
                                SiteModel site,
                                BranchRateModel branch) {
        // sanity check: alignment should have same #taxa as tree
        if (aln.getTaxonSize() != tree.getTree().getLeafCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
        boolean forceJava = Boolean.valueOf(System.getProperty("java.only"));
        if (forceJava) return;
        _alignment = aln;
        _geneTree = tree;
        _nNodeCount = _geneTree.getTree().getNodeCount();
        _siteModel = site;
        _siteModel.setDataType(_alignment.getDataType());
        _branchRateModel = branch;
        if(_branchRateModel == null) {
            _branchRateModel = new StrictClockModel();
        }
        _substitutionModel = _siteModel.getSubstitutionModel();

        _prevBranchLengths = new double[_nNodeCount];
        _storedBranchLengths = new double[_nNodeCount];

        _nStateCount = _alignment.getMaxStateCount();
        _patternCount = _alignment.getPatternCount();
        _eigenCount = 1;

        double[] categoryRates = _siteModel.getCategoryRates(null, null);
        // check for invariant rates category
        if (_siteModel._hasPropInvariantCategory) {
            for (int i = 0; i < categoryRates.length; i++) {
                if (categoryRates[i] == 0) {
                    _proportionInvariant = _siteModel.getRateForCategory(i, _geneTree, null);
                    int stateCount = _alignment.getMaxStateCount();
                    int patterns = _alignment.getPatternCount();
                    calcConstantPatternIndices(patterns, stateCount);
                    _invariantCategory = i;

                    double [] tmp = new double [categoryRates.length - 1];
                    for (int k = 0; k < _invariantCategory; k++) {
                        tmp[k] = categoryRates[k];
                    }
                    for (int k = _invariantCategory + 1; k < categoryRates.length; k++) {
                        tmp[k-1] = categoryRates[k];
                    }
                    categoryRates = tmp;
                    break;
                }
            }
            if (_constantPattern != null && _constantPattern.size() > _alignment.getPatternCount()) {
                // if there are many more constant patterns than patterns (each pattern can
                // have a number of constant patters, one for each state) it is less efficient
                // to just calculate the TreeLikelihood for constant sites than optimising
                System.err.println("switch off constant sites optimisiation: calculating through separate TreeLikelihood category (as in the olden days)");
                _invariantCategory = -1;
                _proportionInvariant = 0;
                _constantPattern = null;
                categoryRates = _siteModel.getCategoryRates(null, null);
            }
        }

        this._categoryCount = _siteModel.getCategoryCount() - (_invariantCategory >= 0 ? 1 : 0);
        _tipCount = _geneTree.getTree().getLeafCount();
        _internalNodeCount = _nNodeCount - _tipCount;

        int compactPartialsCount = _tipCount;
        if (_useAmbiguities) {
            compactPartialsCount = 0;
        }
        // one partials buffer for each tip and two for each internal node (for store restore)
        _partialBufferHelper = new BufferIndexHelper(_nNodeCount, _tipCount);
        // two eigen buffers for each decomposition for store and restore.
        _eigenBufferHelper = new BufferIndexHelper(_eigenCount, 0);
        // two matrices for each node less the root
        _matrixBufferHelper = new BufferIndexHelper(_nNodeCount, 0);
        // one scaling buffer for each internal node plus an extra for the accumulation, then doubled for store/restore
        _scaleBufferHelper = new BufferIndexHelper(getScaleBufferCount(), 0);

        // Attempt to get the resource order from the System Property
        if (_resourceOrder == null) {
            _resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
        }
        if (_preferredOrder == null) {
            _preferredOrder = parseSystemPropertyIntegerArray(PREFERRED_FLAGS_PROPERTY);
        }
        if (_requiredOrder == null) {
            _requiredOrder = parseSystemPropertyIntegerArray(REQUIRED_FLAGS_PROPERTY);
        }
        if (_scalingOrder == null) {
            _scalingOrder = parseSystemPropertyStringArray(SCALING_PROPERTY);
        }

        // first set the rescaling scheme to use from the parser
        _rescalingScheme = PartialsRescalingScheme.DEFAULT;// = rescalingScheme;
        _rescalingScheme = DEFAULT_RESCALING_SCHEME;
        int[] resourceList = null;
        long preferenceFlags = 0;
        long requirementFlags = 0;

        if (_scalingOrder.size() > 0) {
            this._rescalingScheme = PartialsRescalingScheme.parseFromString(
                    _scalingOrder.get(_instanceCount % _scalingOrder.size()));
        }
        if (_resourceOrder.size() > 0) {
            // added the zero on the end so that a CPU is selected if requested resource fails
            resourceList = new int[]{_resourceOrder.get(_instanceCount % _resourceOrder.size()), 0};
            if (resourceList[0] > 0) {
                preferenceFlags |= BeagleFlag.PROCESSOR_GPU.getMask(); // Add preference weight against CPU
            }
        }
        if (_preferredOrder.size() > 0) {
            preferenceFlags = _preferredOrder.get(_instanceCount % _preferredOrder.size());
        }
        if (_requiredOrder.size() > 0) {
            requirementFlags = _requiredOrder.get(_instanceCount % _requiredOrder.size());
        }
        if (_scaling.equals(Scaling.always)) {
            this._rescalingScheme = PartialsRescalingScheme.ALWAYS;
        }
        if (_scaling.equals(Scaling.none)) {
            this._rescalingScheme = PartialsRescalingScheme.NONE;
        }
        if (this._rescalingScheme == PartialsRescalingScheme.DEFAULT) {
            // Define default behaviour here
            //if GPU: the default is^H^Hwas dynamic scaling in BEAST, now NONE
            if (resourceList != null && resourceList[0] > 1) {
                //this.rescalingScheme = PartialsRescalingScheme.DYNAMIC;
                this._rescalingScheme = PartialsRescalingScheme.NONE;
            } else { // if CPU: just run as fast as possible
                //this.rescalingScheme = PartialsRescalingScheme.NONE;
                // Dynamic should run as fast as none until first underflow
                this._rescalingScheme = PartialsRescalingScheme.DYNAMIC;
            }
        }
        if (this._rescalingScheme == PartialsRescalingScheme.AUTO) {
            preferenceFlags |= BeagleFlag.SCALING_AUTO.getMask();
            _useAutoScaling = true;
        }
        String r = System.getProperty(RESCALE_FREQUENCY_PROPERTY);
        if (r != null) {
            _rescalingFrequency = Integer.parseInt(r);
            if (_rescalingFrequency < 1) {
                _rescalingFrequency = RESCALE_FREQUENCY;
            }
        }
        if (preferenceFlags == 0 && resourceList == null) {
            if (_nStateCount == 4 && _patternCount < 10000)
                preferenceFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
        }
        if (_substitutionModel.canReturnComplexDiagonalization()) {
            requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
        }
        _instanceCount++;
        try {
            _beagle = BeagleFactory.loadBeagleInstance(
                    _tipCount,
                    _partialBufferHelper.getBufferCount(),
                    compactPartialsCount,
                    _nStateCount,
                    _patternCount,
                    _eigenBufferHelper.getBufferCount(),            // eigenBufferCount
                    _matrixBufferHelper.getBufferCount(),
                    _categoryCount,
                    _scaleBufferHelper.getBufferCount(), // Always allocate; they may become necessary
                    resourceList,
                    preferenceFlags,
                    requirementFlags
            );
        } catch (Exception e) {
            throw new RuntimeException("Cannot initialize BEAGLE");
        }

        InstanceDetails instanceDetails = _beagle.getDetails();
        ResourceDetails resourceDetails = null;

        if (instanceDetails != null) {
            resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());
            if (resourceDetails != null) {
                StringBuilder sb = new StringBuilder("  Using BEAGLE version: " + BeagleInfo.getVersion()
                        + " resource ");
                sb.append(resourceDetails.getNumber()).append(": ");
                sb.append(resourceDetails.getName()).append("\n");
                if (resourceDetails.getDescription() != null) {
                    String[] description = resourceDetails.getDescription().split("\\|");
                    for (String desc : description) {
                        if (desc.trim().length() > 0) {
                            sb.append("    ").append(desc.trim()).append("\n");
                        }
                    }
                }
                sb.append("    with instance flags: ").append(instanceDetails.toString());
                if (_print) System.out.println(sb.toString());
            } else {
                System.err.println("  Error retrieving BEAGLE resource for instance: " + instanceDetails.toString());
                throw new RuntimeException("Cannot initialize BEAGLE");
            }
        } else {
            System.err.println("  No external BEAGLE resources available, or resource list/requirements not met, using Java implementation");
            throw new RuntimeException("Cannot initialize BEAGLE");
        }

        TNode[] nodes = _geneTree.getNodeArray();
        for (int i = 0; i < _tipCount; i++) {
            int taxon = _alignment.getTaxonIndex(nodes[i].getName());
            if (_useAmbiguities || _useTipLikelihoods) {
                setPartials(_beagle, i, taxon);
            } else {
                setStates(_beagle, i, taxon);
            }
        }

        if (_alignment._isAscertained) {
            _useAscertainedSitePatterns = true;
        }

        double[] patternWeights = new double[_patternCount];
        for (int i = 0; i < _patternCount; i++) {
            patternWeights[i] = _alignment.getPatternWeight(i);
        }
        _beagle.setPatternWeights(patternWeights);

        if (this._rescalingScheme == PartialsRescalingScheme.AUTO &&
                resourceDetails != null &&
                (resourceDetails.getFlags() & BeagleFlag.SCALING_AUTO.getMask()) == 0) {
            // If auto scaling in BEAGLE is not supported then do it here
            this._rescalingScheme = PartialsRescalingScheme.DYNAMIC;
            if (_print) {
                System.out.println("  Auto rescaling not supported in BEAGLE, using : " + this._rescalingScheme.getText());
            }
        } else if (_print) {
            System.out.println("  Using rescaling scheme : " + this._rescalingScheme.getText());
        }

        if (this._rescalingScheme == PartialsRescalingScheme.DYNAMIC) {
            _everUnderflowed = false; // If false, BEAST does not rescale until first under-/over-flow.
        }

        _updateSubstitutionModel = true;
        _updateSiteModel = true;
        setUpSubstModel();
        _beagle.setCategoryRates(categoryRates);
        _currentCategoryRates = categoryRates;
        _currentFreqs = new double[_nStateCount];
        _currentCategoryWeights = new double[categoryRates.length];
    }

    public int getPatternCount() {
        return _patternCount;
    }

    public void setUpSubstModel() {
        for (int i = 0; i < _eigenCount; i++) {
            EigenDecomposition ed = _substitutionModel.getEigenDecomposition();
            _eigenBufferHelper.flipOffset(i);
            _beagle.setEigenDecomposition(
                    _eigenBufferHelper.getOffsetIndex(i),
                    ed.getEigenVectors(),
                    ed.getInverseEigenVectors(),
                    ed.getEigenValues());
        }
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogP() {
        if (_patternLogLikelihoods == null) {
            _patternLogLikelihoods = new double[_patternCount];
        }
        if (_matrixUpdateIndices == null) {
            _matrixUpdateIndices = new int[_eigenCount][_nNodeCount];
            _branchLengths = new double[_eigenCount][_nNodeCount];
            _branchUpdateCount = new int[_eigenCount];
            _scaleBufferIndices = new int[_internalNodeCount];
            _storedScaleBufferIndices = new int[_internalNodeCount];
        }
        if (_operations == null) {
            _operations = new int[1][_internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
            _operationCount = new int[1];
        }
        _recomputeScaleFactors = false;
        if (this._rescalingScheme == PartialsRescalingScheme.ALWAYS) {
            _useScaleFactors = true;
            _recomputeScaleFactors = true;
        } else if (this._rescalingScheme == PartialsRescalingScheme.DYNAMIC && _everUnderflowed) {
            _useScaleFactors = true;
            if (_rescalingCountInner < RESCALE_TIMES) {
                _recomputeScaleFactors = true;
                _hasDirt = 2;
            }
            _rescalingCountInner++;
            _rescalingCount++;
            if (_rescalingCount > RESCALE_FREQUENCY) {
                _rescalingCount = 0;
                _rescalingCountInner = 0;
            }
        } else if (this._rescalingScheme == PartialsRescalingScheme.DELAYED && _everUnderflowed) {
            _useScaleFactors = true;
            _recomputeScaleFactors = true;
            _hasDirt = 2;
            _rescalingCount++;
        }
        for (int i = 0; i < _eigenCount; i++) {
            _branchUpdateCount[i] = 0;
        }
        _operationListCount = 0;
        _operationCount[0] = 0;
        final TNode root = _geneTree.getTree().getRoot();
        traverse(root, null, true);

        if (_updateSubstitutionModel) {
            setUpSubstModel();
        }
        if (_updateSiteModel) {
            double[] categoryRates = _siteModel.getCategoryRates(null, null);
            if (_constantPattern != null) {
                double [] tmp = new double [categoryRates.length - 1];
                for (int k = 0; k < _invariantCategory; k++) {
                    tmp[k] = categoryRates[k];
                }
                for (int k = _invariantCategory + 1; k < categoryRates.length; k++) {
                    tmp[k-1] = categoryRates[k];
                }
                categoryRates = tmp;
            }
            for (int i = 0; i < categoryRates.length; i++) {
                if (categoryRates[i] != _currentCategoryRates[i]) {
                    _beagle.setCategoryRates(categoryRates);
                    i = categoryRates.length;
                }
            }
            _currentCategoryRates = categoryRates;
        }

        for (int i = 0; i < _eigenCount; i++) {
            if (_branchUpdateCount[i] > 0) {
                _beagle.updateTransitionMatrices(
                        _eigenBufferHelper.getOffsetIndex(i),
                        _matrixUpdateIndices[i],
                        null,
                        null,
                        _branchLengths[i],
                        _branchUpdateCount[i]);
            }
        }
        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;
        do {
            _beagle.updatePartials(_operations[0], _operationCount[0], Beagle.NONE);
            int rootIndex = _partialBufferHelper.getOffsetIndex(_geneTree.getNodeLabel(root));
            double[] categoryWeights = _siteModel.getCategoryProportions(null, null);
            if (_constantPattern != null) {
                double [] tmp = new double [categoryWeights.length - 1];
                for (int k = 0; k < _invariantCategory; k++) {
                    tmp[k] = categoryWeights[k];
                }
                for (int k = _invariantCategory + 1; k < categoryWeights.length; k++) {
                    tmp[k-1] = categoryWeights[k];
                }
                categoryWeights = tmp;
            }
            double[] frequencies = _substitutionModel.getFrequencies();

            int cumulateScaleBufferIndex = Beagle.NONE;
            if (_useScaleFactors) {
                if (_recomputeScaleFactors) {
                    _scaleBufferHelper.flipOffset(_internalNodeCount);
                    cumulateScaleBufferIndex = _scaleBufferHelper.getOffsetIndex(_internalNodeCount);
                    _beagle.resetScaleFactors(cumulateScaleBufferIndex);
                    _beagle.accumulateScaleFactors(_scaleBufferIndices, _internalNodeCount, cumulateScaleBufferIndex);
                } else {
                    cumulateScaleBufferIndex = _scaleBufferHelper.getOffsetIndex(_internalNodeCount);
                }
            } else if (_useAutoScaling) {
                _beagle.accumulateScaleFactors(_scaleBufferIndices, _internalNodeCount, Beagle.NONE);
            }

            // these could be set only when they change but store/restore would need to be considered
            for (int i = 0; i < categoryWeights.length; i++) {
                if (categoryWeights[i] != _currentCategoryWeights[i]) {
                    _beagle.setCategoryWeights(0, categoryWeights);
                    i = categoryWeights.length;
                }
            }
            _currentCategoryWeights = categoryWeights;
            for (int i = 0; i < frequencies.length; i++) {
                if (frequencies[i] != _currentFreqs[i]) {
                    _beagle.setStateFrequencies(0, frequencies);
                    i = frequencies.length;
                }
            }
            _currentFreqs = frequencies;

            double[] sumLogLikelihoods = new double[1];
            _beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                    new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

            logL = sumLogLikelihoods[0];

            if (_useAscertainedSitePatterns) {
                _beagle.getSiteLogLikelihoods(_patternLogLikelihoods);
                logL = getAscertainmentCorrectedLogLikelihood(_alignment,
                        _patternLogLikelihoods, _alignment.getPatternWeights(), frequencies);
            } else if (_invariantCategory >= 0) {
                _beagle.getSiteLogLikelihoods(_patternLogLikelihoods);
                int [] patternWeights = _alignment.getPatternWeights();
                _proportionInvariant = _siteModel.getProportionInvariant();

                for (int k : _constantPattern) {
                    int i = k / _nStateCount;
                    int j = k % _nStateCount;
                    _patternLogLikelihoods[i] = (Math.log(Math.exp(_patternLogLikelihoods[i])
                            + _proportionInvariant * frequencies[j]));
                }
                logL = 0.0;
                for (int i = 0; i < _patternCount; i++) {
                    logL += _patternLogLikelihoods[i] * patternWeights[i];
                }
            }

            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                _everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt && (_rescalingScheme == PartialsRescalingScheme.DYNAMIC
                        || _rescalingScheme == PartialsRescalingScheme.DELAYED)) {
                    _useScaleFactors = true;
                    _recomputeScaleFactors = true;
                    for (int i = 0; i < _eigenCount; i++) {
                        _branchUpdateCount[i] = 0;
                    }
                    _operationCount[0] = 0;
                    traverse(root, null, false);
                    done = false;
                    firstRescaleAttempt = false;
                } else {
                    done = true;
                }
            } else {
                done = true;
            }

        } while (!done);

        _updateSubstitutionModel = false;
        _updateSiteModel = false;
        _logP = logL;
        return logL;
    }

    public List<String> getArguments() {return null;}

    public List<StateNode> getConditions() {
        return _siteModel.getConditions();
    }


    @Override
    public double propose() {
        return 0;
    }

    @Override
    public void undo() {

    }

    @Override
    public void accept() {

    }

    @Override
    public void reject() {

    }

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

    // --- private/protected methods ---

    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @param node           node
     * @param operatorNumber operatorNumber
     * @param flip           flip
     * @return boolean
     */
    private int traverse(TNode node, int[] operatorNumber, boolean flip) {

        int nodeNum = _geneTree.getNodeLabel(node);
        if (operatorNumber != null) {
            operatorNumber[0] = -1;
        }
        // First update the transition probability matrix(ices) for this branch
        int update = _hasDirt;
        final double branchRate = _branchRateModel.getRateForBranch(_geneTree, node);
        final double branchTime = node.getParentDistance() * branchRate;
        if (!node.isRoot() && (update != 0 || branchTime != _prevBranchLengths[nodeNum])) {
            _prevBranchLengths[nodeNum] = branchTime;
            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }
            if (flip) {
                // first flip the matrixBufferHelper
                _matrixBufferHelper.flipOffset(nodeNum);
            }
            // then set which matrix to update
            final int eigenIndex = 0;// = m_substitutionModel.getBranchIndex(node);
            final int updateCount = _branchUpdateCount[eigenIndex];
            _matrixUpdateIndices[eigenIndex][updateCount] = _matrixBufferHelper.getOffsetIndex(nodeNum);
            _branchLengths[eigenIndex][updateCount] = branchTime;
            _branchUpdateCount[eigenIndex]++;
            update |= 1;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {
            // Traverse down the two child nodes
            List<TNode> children = IterableHelp.toList(node.getChildren());

            TNode child1 = children.get(0);
            final int[] op1 = {-1};
            final int update1 = traverse(child1, op1, flip);

            TNode child2 = children.get(1);
            final int[] op2 = {-1};
            final int update2 = traverse(child2, op2, flip);

            // If either child node was updated then update this node too
            if (update1 != 0 || update2 != 0) {
                int x = _operationCount[_operationListCount] * Beagle.OPERATION_TUPLE_SIZE;
                if (flip) {
                    // first flip the partialBufferHelper
                    _partialBufferHelper.flipOffset(nodeNum);
                }
                final int[] operations = _operations[_operationListCount];
                operations[x] = _partialBufferHelper.getOffsetIndex(nodeNum);
                if (_useScaleFactors) {
                    // get the index of this scaling buffer
                    int n = nodeNum - _tipCount;
                    if (_recomputeScaleFactors) {
                        // flip the indicator: can take either n or (internalNodeCount + 1) - n
                        _scaleBufferHelper.flipOffset(n);
                        // store the index
                        _scaleBufferIndices[n] = _scaleBufferHelper.getOffsetIndex(n);
                        operations[x + 1] = _scaleBufferIndices[n]; // Write new scaleFactor
                        operations[x + 2] = Beagle.NONE;
                    } else {
                        operations[x + 1] = Beagle.NONE;
                        operations[x + 2] = _scaleBufferIndices[n]; // Read existing scaleFactor
                    }
                } else {
                    if (_useAutoScaling) {
                        _scaleBufferIndices[nodeNum - _tipCount] = _partialBufferHelper.getOffsetIndex(nodeNum);
                    }
                    operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                    operations[x + 2] = Beagle.NONE;
                }
                operations[x + 3] = _partialBufferHelper.getOffsetIndex(_geneTree.getNodeLabel(child1)); // source node 1
                operations[x + 4] = _matrixBufferHelper.getOffsetIndex(_geneTree.getNodeLabel(child1)); // source matrix 1
                operations[x + 5] = _partialBufferHelper.getOffsetIndex(_geneTree.getNodeLabel(child2)); // source node 2
                operations[x + 6] = _matrixBufferHelper.getOffsetIndex(_geneTree.getNodeLabel(child2)); // source matrix 2
                _operationCount[_operationListCount]++;
                update |= (update1 | update2);
            }
        }
        return update;
    }

    private void calcConstantPatternIndices(final int patterns, final int stateCount) {
        _constantPattern = new ArrayList<>();
        for (int i = 0; i < patterns; i++) {
            final int[] pattern = _alignment.getPattern(i);
            final boolean[] isInvariant = new boolean[stateCount];
            Arrays.fill(isInvariant, true);
            for (final int state : pattern) {
                final boolean[] isStateSet = _alignment.getStateSet(state);
                if (_useAmbiguities || !_alignment.getDataType().isAmbiguousState(state)) {
                    for (int k = 0; k < stateCount; k++) {
                        isInvariant[k] &= isStateSet[k];
                    }
                }
            }
            for (int k = 0; k < stateCount; k++) {
                if (isInvariant[k]) {
                    _constantPattern.add(i * stateCount + k);
                }
            }
        }
    }

    /**
     * Sets the partials from a sequence in an alignment.
     */
    protected final void setPartials(Beagle beagle,
                                     int nodeIndex, int taxon) {

        double[] partials = new double[_patternCount * _nStateCount * _categoryCount];

        int v = 0;
        for (int i = 0; i < _patternCount; i++) {
            double[] tipProbabilities = _alignment.getTipLikelihoods(taxon, i);
            if (tipProbabilities != null) {
                for (int state = 0; state < _nStateCount; state++) {
                    partials[v++] = tipProbabilities[state];
                }
            } else {
                int stateCount = _alignment.getPattern(taxon, i);
                boolean[] stateSet = _alignment.getStateSet(stateCount);
                for (int state = 0; state < _nStateCount; state++) {
                    partials[v++] = (stateSet[state] ? 1.0 : 0.0);
                }
            }
        }
        int n = _patternCount * _nStateCount;
        int k = n;
        for (int i = 1; i < _categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }
        beagle.setPartials(nodeIndex, partials);
    }

    /**
     * Sets the partials from a sequence in an alignment.
     */
    protected final void setStates(Beagle beagle,
                                   int nodeIndex, int taxon) {

        int[] states = new int[_patternCount];
        for (int i = 0; i < _patternCount; i++) {
            int state = _alignment.getPattern(taxon, i);
            states[i] = state;
        }
        beagle.setTipStates(nodeIndex, states);
    }

    protected void setPartials(int number, double[] partials) {
        _beagle.setPartials(_partialBufferHelper.getOffsetIndex(number), partials);
    }

    private double getAscertainmentCorrectedLogLikelihood(Alignment patternList,
                                                          double[] patternLogLikelihoods,
                                                          int[] patternWeights,
                                                          double [] frequencies) {
        if (_constantPattern != null) {
            _proportionInvariant = _siteModel.getProportionInvariant();
            for (int k : _constantPattern) {
                int i = k / _nStateCount;
                int j = k % _nStateCount;
                patternLogLikelihoods[i] = (Math.log(Math.exp(patternLogLikelihoods[i])
                        + _proportionInvariant * frequencies[j]));
            }
        }
        double logL = 0.0;
        double ascertainmentCorrection = patternList.getAscertainmentCorrection(patternLogLikelihoods);
        for (int i = 0; i < _patternCount; i++) {
            logL += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
        }
        return logL;
    }

    private int getScaleBufferCount() {
        return _internalNodeCount + 1;
    }

    private static List<Integer> parseSystemPropertyIntegerArray(String propertyName) {
        List<Integer> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    private static List<String> parseSystemPropertyStringArray(String propertyName) {

        List<String> order = new ArrayList<>();

        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    String s = part.trim();
                    order.add(s);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid getEigenDecompositionentry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    // test
    public static void main(String[] args) {
        String newick1 = "((C:0.002859330148428049,G:0.002859330148428049):0.010473421735219215,(R:0.010936563366414687,(L:0.008340943025951807,(A:0.004491653225077681,Q:0.004491653225077681):0.0038492898008741254):0.00259562034046288):0.0023961885172325784);";
        // -2125.1658868743093
        String newick2 = "((L:0.087137,R:0.087137):0.036282,((G:0.037923,C:0.037923):0.069297,(Q:0.032283,A:0.032283):0.074936):0.016199):0;";
        // -1731.9586701487046
        String newick3 = "((C:0.002859330148428049,G:0.002859330148428049)I3:0.11047342173521922,(R:0.010936563366414687,(L:0.008340943025951807,(A:0.004491653225077681,Q:0.004491653225077681)I0:0.0038492898008741254)I1:0.00259562034046288)I2:0.10239618851723259)I4;";
        // -2085.91597350229
        Map<String, String> omap = new HashMap<>();
        omap.put("L", "GTCATCCGACTCGTCTGACCCCGGCGCGGATCTGGAACCCCTGCCGTGACTCGGCGCGCAGGCGGGCGCATCGCCAAACCTAGCCCAGTGCAGGTCGTGGATGTAACGCTAGCGCTATTTGTGCATCAAAAGATAGCCGCTCGACGTCGCGCCCTTGTGGCCTTACCATCAGTGCGGCCCATCGGCCAGACCTCTGGAAGACTGAATGACAATATTTCGAAGGCGGAAGCCCGAAATGGCGTATGTGTCTGGAACAGAGGAGGGGTAGACTGCCCGTATGCACAGAACGCGACTGTGCGCTAGGCAAGGCTTGACCTAGCCGAATCGGTGCAGACTTTATTGCGACCGCAACGCTTAAGCCCCAGTGCGTCATGCCCCGTAAGTGTATGGGGCACTGCCCTTCGGCGCCCATCGTTGCCGCCTGGACGGCGACTCTCCTCATGGATTCCGCGCAACGGCCGGCTACCGCGAACGAGCCGTCATACACGTTTCGTACAGTA");
        omap.put("R", "GTTCACCGACTCGTCTGATCCCGACGCGGATCTGGGACCCATGCCATGGCTCAAGGCGCAGGCGGACGCATCGCCAAACCTAGCCCAGTGTAGGTCGTGGATGTCACACTAGGGCTATTTGTGCATCAGAAGATAGCTGCTCGACCTAGCGCCCTTGTTGCCCTGCCATCGGTGCGACCCATCGACCAGATCTCTGGAAGGCTGAACGATCATCTTGCGAAGGCGGAAGCCCGAAGTAACGTAGGTATCTGTAGTAGAAAAAGGGTAGACTACCCGTGTGCATAGAGTGCGACTGCGCGCTAGGCAAAGCTTGAACTAGCCGAATCGGTGCAGACCTTGTTGCGACCGCAATGCTTATGCCCCTGAGCGTCATGCGTCGTAAGTGTGTAAGGCACTGCGCTTCGGCGCCCCTCGCTGCCGCCTGGACGGCGCCTCTCCTCATGAATTCCGCGCAACGGCCGACTACCGCGAACGAGCCGCCATACACTTTTCGCGCATTG");
        omap.put("G", "GTCCTCCGACTCGTCTGATCCCGACGCGAACCCGGCGCGCCGGCCATGACTCGGGGCGCAGGCGGGCGCACCGTCAAGCCTAGCCCAGGGCAGATCGTGGATATCACGCTAGCGCTATTTGTGTATCAAAGCATAACTGCTCGGCGTAGCACCCTCGTGGCCTCTTCATCAGCGCGGCCCGTCGACCAGATTTCTGGAAGGCTGGATGACCATTTTGCGAGGGCGGGAGTCCGAAGTAGCGTAAGTGTCCGGAGCAGAGAAGGAGTAAACTGCCGGTGTGCTTAGAGCGCGACTGTGCGCCAGGCAAAACTTGAAATAGCCGAATCGATGCACACTCCGTTGCGACCGCAGTGCTTATGCCCCTGTGCGTCCTGCATCGTAAGGGTATAGGGCGCGGCGCTTCAGCGCCCGTCGTTGCCGCCTGGACGACGCCCCTTCTCACGGATGCCGCCCAACGGCCCGCTACCCCGAACGGTCAGCCAGACACTTTTCGTACATTA");
        omap.put("C", "GCCCTCCGACTCGTCTGATCCCGACGCGGACCCGACGCACCGGCCATGACTCGGGGCGCAGGCGGGCGCACCGTCAAGCCTAGCCCAGGGCAGACCGTGGATATCGCGCTAGCGCTATTTGTGCATCAAAACGTAACTGCTCGGCGTAGCACCCTTGTGGCCGCTCCGTCAGCGCGGCCCGTCGACCAGATTTCTGGAAGGCTGGATGACCATCTTGCGGGGGCGGGAGTCCGCAGTAGCGTAAGTGTCCGGAGCAGAGAAGGAGTAAACTGCCGGTATTCATAAAGTGCAACTGTGCGCCAGGCAAAACTTGAAATAGCCGATTCGATGCACACCTCGTTGCAACCGAAGTGCCTATGCCCCTGTGCGTCCTGCATCGTAAGGGTATAGGGCGCGGCGCTTCAGCGCCCACCGTTGCCGTCTGAACGGCGCCTCTTCTCACGGATGCCGCCCAACGGCCGGCTACCCCGAACGGTCAGCCAGGCACTTTTCGTACATTA");
        omap.put("Q", "GTCCTTCGACTCGTCTGATCCCGACGCGGACCCGGCGCCCCTGCCATGACTCGGGGCGCAGGCAGGCGCACCGTCAAATCGAGCCCAGTGCAGGTCGTGGATGTCACGCTAGCGCTATCTGTGCATCAAGGGGTAACTGTTCGCCGTAGTGCCCTTGTGGCCTCGCCATCAGTGCGGCCCATCGACCAGGTCTCTGGAAGGCTGGATAACCACTTTGCGAAGACGGGAGTCCGAAGTGGCGTATGTGTCCAGAGCAGAGAAGGGATAAACTGCCCGAATGCGTAGAGCGCGACCGTGCGCTAGGCAAAGTCTGAGCTAGCCGAATCGGCGCGCACTTTATTGCGGCCGCAGTGCCTTGGTCCCTGTGCGTCCCGCAACGTAAGGGTATGGGGCACCGCGCTTCGGCGCCTATCCTTGCCGCCTGGACGGCGCTTCTTCTCATGGATTCCACGCAACGACCGGCACCTGCGAACGGTCCGCCAGACACTTTTCGTACACTA");
        omap.put("A", "GTCCTCCGACTCATCTGATCCCGACGCGGACCCGGCGCCCCGGCCATGACTCGGGGCGCAGGCGGGCGCACCGTCAAGCCTTGCCCGGTGCAGGTCGTGGATGTCACGCTAGCGCTATCTGCGCATCAAAAGGTAACTGTTCGACGTGGTGCCCTTGTGGCCTCGCCATCAGTGCGGCCCATCGACCAGGTCTCTGGAAGGCTGTATAACCACTTTGCGACGACGGGAGTCCGAAGTGGCGTAAGTGTCCGGAGCAGAGAAGGGATAAACTGCCTGAATGCGCAGAGCGCGACCGTGCGCTAGGCAAAGTTTGAGCTAGCCGAATCGGCGCGCACTTTATTGCGGCCGCAGTGTGTTGGTCCCTGTGCGTCCCGCAACGTAAGGGTATGGGGCACCGCGCTTCGGCGCCTATCGTTGCCGCCTGTACGGCGATTCTTCTCATGAATCCCACACATCGACCGGCACCTGCGAACGGTCCGCCAGACACTTTTCGTACACTA");

        Alignment aln = new Alignment(omap);

        double[] base = {0.2112, 0.2888, 0.2896, 0.2104};
        double[] trans = {0.2173/0.2070, 0.9798/0.2070, 0.2575/0.2070, 0.1038/0.2070, 1/0.2070, 1.0};
        Utils._BASE_FREQS = base;
        Utils._TRANS_RATES = base;
        Utils._SUBSTITUTION_MODEL = "GTR";

        UltrametricTree uTree1 = new UltrametricTree(newick1, aln);
        Trees.autoLabelNodes((MutableTree) uTree1.getTree());
        System.out.println(uTree1.toString());
        for(int i = 0; i < 1000; i++) {
            uTree1.propose();
        }
        System.out.println(uTree1.toString());
        newick3 = uTree1.toString();
        UltrametricTree uTree2 = new UltrametricTree(newick3, aln);
        System.out.println(uTree2.toString());

        // use aln, tree, site = new SiteModel(new GTR()), null
        Frequencies freq = new Frequencies(aln, base);
        SubstitutionModel subst = new GTR(freq, trans);

        BeagleTreeLikelihood beagle1 = new BeagleTreeLikelihood(aln, uTree1, new SiteModel(subst), null);
        System.out.println(beagle1._hasDirt);
        System.out.println(beagle1.calculateLogP());
        System.out.println(beagle1._hasDirt);
        TNode root = uTree1.getTree().getRoot();
        uTree1.setNodeHeight(root, uTree1.getNodeHeight(root) + 0.1);
        System.out.println(uTree1);
        System.out.println(beagle1.calculateLogP());

        BeagleTreeLikelihood beagle2 = new BeagleTreeLikelihood(aln, uTree2, new SiteModel(subst), null);
        System.out.println(beagle2._hasDirt);
        System.out.println(beagle2.calculateLogP());
        System.out.println(beagle2._hasDirt);
    }
}
