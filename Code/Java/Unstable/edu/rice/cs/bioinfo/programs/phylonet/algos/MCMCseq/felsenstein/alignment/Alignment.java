package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.datatype.DataType;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.datatype.Nucleotide;

import java.util.*;

/**
 * Created by wendingqiao on 5/3/16.
 */
public class Alignment {

    protected List<String> _taxaNames = new ArrayList<>(); // sorted
    protected List<Sequence> _sequences = new ArrayList<>(); // corresponds to _taxaNames

    // state counts for each sequence
    protected List<Integer> _stateCounts = new ArrayList<>();
    protected int _maxStateCount = -1;
    // state codes for the sequences, a matrix of #taxa X #sites
    protected List<List<Integer>> _counts = new ArrayList<>();

    protected DataType m_dataType = new Nucleotide();

    // weights of sites, default value is 1
    protected int[] _siteWeights = null;
    // weight over the columns of a matrix
    protected int[] _patternWeight;

    // Probabilities associated with leaves when the characters are uncertain
    public List<double[][]> _tipLikelihoods = new ArrayList<>(); // #taxa X #sites X #states
    private boolean _usingTipLikelihoods = false;

    // pattern state encodings, a matrix of #patterns X #taxa
    protected int [][] _sitePatterns;
    // maps site nr to pattern nr, an array of #sites
    protected int[] _patternIndex;

    // From AscertainedAlignment
    public boolean _isAscertained = false;
    protected Set<Integer> _excludedPatterns;
    protected int _excludefrom = 0; // first site to condition on
    protected int _excludeto = 0; // last site to condition on (but excluding this site
    protected int _excludeevery = 1; // interval between sites to condition on

    private Map<String, String> _aln;
    private String _name = null;

    public Alignment(Map<String, String> sequences, String name) {
        this(sequences);
        this._name = name;
    }

    public Alignment(Map<String, String> sequences) {
        _aln = sequences;
        for(String key : sequences.keySet()) {
            if (_taxaNames.contains(key)) {
                throw new RuntimeException("Duplicate taxon found in alignment: " + key);
            }
            _taxaNames.add(key);
        }
        Collections.sort(_taxaNames);

        for(String taxon: _taxaNames) {
            Sequence seq = new Sequence(taxon, sequences.get(taxon));
            _sequences.add(seq);
            _counts.add(seq.getSequence(m_dataType));
            if(seq.getStateCount() == -1) {
                throw new RuntimeException("state count has not been initialized yet " + taxon + " " + sequences.get(taxon));
            }
            _stateCounts.add(seq.getStateCount());
            _maxStateCount = Math.max(_maxStateCount, seq.getStateCount());
        }
        if (_counts.size() == 0) {
            throw new RuntimeException("Sequence data expected, but none found");
        }
        calcPatterns();
    }

    public Map<String, String> getAlignment() {
        return _aln;
    }

    public String getName() { return _name; }

    private void setupAscertainment(int from, int to, int every) {
        _isAscertained = true;
        _excludefrom = from;
        _excludeto = to;
        _excludeevery = every;
        _excludedPatterns = new HashSet<>();
        for (int i = from; i < to; i += every) {
            int patternIndex = _patternIndex[i];
            // reduce weight, so it does not confuse the tree likelihood
            _patternWeight[patternIndex] = 0;
            _excludedPatterns.add(patternIndex);
        }
    }


    /**
     * SiteComparator is used for ordering the sites,
     * which makes it easy to identify patterns.
     */
    class SiteComparator implements Comparator<int[]> {
        @Override
        public int compare(int[] o1, int[] o2) {
            for (int i = 0; i < o1.length; i++) {
                if (o1[i] > o2[i]) return 1;
                if (o1[i] < o2[i]) return -1;
            }
            return 0;
        }
    }

    /**
     * calculate patterns from sequence data
     */
    private void calcPatterns() {
        int taxonCount = _counts.size();
        int siteCount = _counts.get(0).size();
        // convert data to transposed int array
        int[][] data = new int[siteCount][taxonCount];
        for (int i = 0; i < taxonCount; i++) {
            List<Integer> sites = _counts.get(i);
            for (int j = 0; j < siteCount; j++) {
                data[j][i] = sites.get(j);
            }
        }
        // sort data
        SiteComparator comparator = new SiteComparator();
        Arrays.sort(data, comparator);
        // count patterns in sorted data
        int patterns = 1;
        int[] weights = new int[siteCount];
        weights[0] = 1;
        for (int i = 1; i < siteCount; i++) {
            if (_usingTipLikelihoods || comparator.compare(data[i - 1], data[i]) != 0) {
                // In the case where we're using tip probabilities, we need to treat each
                // site as a unique pattern, because it could have a unique probability vector.
                patterns++;
                data[patterns - 1] = data[i];
            }
            weights[patterns - 1]++;
        }
        // reserve memory for patterns
        _patternWeight = new int[patterns];
        _sitePatterns = new int[patterns][taxonCount];
        for (int i = 0; i < patterns; i++) {
            _patternWeight[i] = weights[i];
            _sitePatterns[i] = data[i];
        }
        // find patterns for the sites
        _patternIndex = new int[siteCount];
        for (int i = 0; i < siteCount; i++) {
            int[] sites = new int[taxonCount];
            for (int j = 0; j < taxonCount; j++) {
                sites[j] = _counts.get(j).get(i);
            }
            _patternIndex[i] = Arrays.binarySearch(_sitePatterns, sites, comparator);
        }
    }

    public List<String> getTaxaNames() {
        return _taxaNames;
    }

    public List<Integer> getStateCounts() {
        return _stateCounts;
    }

    public DataType getDataType() {
        return m_dataType;
    }

    public List<List<Integer>> getCounts() {
        return _counts;
    }

    public int getTaxonSize() {
        return _taxaNames.size();
    }

    public int getTaxonIndex(String taxon) {
        return _taxaNames.indexOf(taxon);
    }

    public int getMaxStateCount() {
        return _maxStateCount;
    }

    public boolean[] getStateSet(int state) {
        return m_dataType.getStateSet(state);
    }

    public int getPatternCount() {
        return _sitePatterns.length;
    }

    public int[] getPattern(int patternIndex) {
        return _sitePatterns[patternIndex];
    }

    public int getPattern(int taxonIndex, int patternIndex) {
        return _sitePatterns[patternIndex][taxonIndex];
    }

    public int getPatternWeight(int patternIndex) {
        return _patternWeight[patternIndex];
    }

    public int getPatternIndex(int site) {
        return _patternIndex[site];
    }

    public int getSiteCount() {
        return _patternIndex.length;
    }

    public int[] getPatternWeights() {
        return _patternWeight;
    }

    public double[] getTipLikelihoods(int taxon, int idx) { return null; }

    private long getTotalWeight() {
        long totalWeight = 0;
        for (int weight : _patternWeight) {
            totalWeight += weight;
        }
        return totalWeight;
    }


    //Methods from AscertainedAlignment
    public Set<Integer> getExcludedPatternIndices() {
        return _excludedPatterns;
    }

    public int getExcludedPatternCount() {
        return _excludedPatterns.size();
    }

    public double getAscertainmentCorrection(double[] patternLogProbs) {
        double excludeProb = 0, includeProb = 0, returnProb = 1.0;
        for (int i : _excludedPatterns) {
            excludeProb += Math.exp(patternLogProbs[i]);
        }
        if (includeProb == 0.0) {
            returnProb -= excludeProb;
        } else if (excludeProb == 0.0) {
            returnProb = includeProb;
        } else {
            returnProb = includeProb - excludeProb;
        }
        return Math.log(returnProb);
    }


    public static String getSequence(Alignment data, int taxonIndex) {

        int[] states = new int[data.getPatternCount()];
        for (int i = 0; i < data.getPatternCount(); i++) {
            int[] sitePattern = data.getPattern(i);
            states[i] = sitePattern[taxonIndex];
        }
        try {
            return data.getDataType().stateToString(states);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return null;
    }

    // test
    public static void main(String[] args) {
        Map<String, String> input = new HashMap<>();
        input.put("A", "ATATCG");
        input.put("B", "ATATG-");
        input.put("C", "CTAT-G");
        Alignment aln = new Alignment(input);
        System.out.println(aln.getTaxonSize() == 3);
        List<String> taxa = aln.getTaxaNames();
        System.out.println(taxa.get(0) == "A" && taxa.get(1) == "B" && taxa.get(2) == "C");
        System.out.println(aln.getMaxStateCount() == 4);
        System.out.println(aln.getPatternCount() == 5);
        System.out.println(Arrays.toString(aln.getPatternWeights())); // 1,1,1,1,2
        System.out.println(Arrays.toString(aln.getStateSet(0))); // TFFF
        System.out.println(Arrays.toString(aln.getStateSet(17))); // TTTT
        System.out.println(Arrays.toString(aln.getPattern(4))); // T T T -> 3, 3, 3
    }

}
