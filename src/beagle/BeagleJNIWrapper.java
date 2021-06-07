/*
 * BeagleJNIjava
 *
 */

package beagle;


/*
 * BeagleJNIjava
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 *
 */

public class BeagleJNIWrapper {

    public static final String LIBRARY_NAME = getPlatformSpecificLibraryName();

    /**
     * private constructor to enforce singleton instance
     */
    private BeagleJNIWrapper() {
    }

    public native String getVersion();

    public native String getCitation();

    public native ResourceDetails[] getResourceList();

    public native BenchmarkedResourceDetails[] getBenchmarkedResourceList(
                                                    int tipCount,
                                                    int compactBufferCount,
                                                    int stateCount,
                                                    int patternCount,
                                                    int categoryCount,
                                                    final int[] resourceList,
                                                    int resourceCount,
                                                    long preferenceFlags,
                                                    long requirementFlags,
                                                    int eigenModelCount,
                                                    int partitionCount,
                                                    int calculateDerivatives,
                                                    long benchmarkFlags);

    public native int createInstance(
            int tipCount,
            int partialsBufferCount,
            int compactBufferCount,
            int stateCount,
            int patternCount,
            int eigenBufferCount,
            int matrixBufferCount,
            int categoryCount,
            int scaleBufferCount,
            final int[] resourceList,
            int resourceCount,
            long preferenceFlags,
            long requirementFlags,
            InstanceDetails returnInfo);

    public native int finalize(int instance);

    public native int setCPUThreadCount(int instance, int threadCount);

    public native int setPatternWeights(int instance, final double[] patternWeights);

    public native int setPatternPartitions(int instance, int partitionCount, final int[] patternPartitions);

    public native int setTipStates(int instance, int tipIndex, final int[] inStates);

    public native int getTipStates(int instance, int tipIndex, final int[] inStates);

    public native int setTipPartials(int instance, int tipIndex, final double[] inPartials);

    public native int setPartials(int instance, int bufferIndex, final double[] inPartials);

    public native int getPartials(int instance, int bufferIndex, int scaleIndex,
                                  final double[] outPartials);
    
    public native int getLogScaleFactors(int stance, int scaleIndex, final double[] outFactors);

    public native int setEigenDecomposition(int instance,
                                            int eigenIndex,
                                            final double[] eigenVectors,
                                            final double[] inverseEigenValues,
                                            final double[] eigenValues);

    public native int setStateFrequencies(int instance,
                                          int stateFrequenciesIndex,
                                          final double[] stateFrequencies);

    public native int setCategoryWeights(int instance,
                                         int categoryWeightsIndex,
                                         final double[] categoryWeights);

    public native int setCategoryRates(int instance,
                                       final double[] inCategoryRates);

    public native int setCategoryRatesWithIndex(int instance,
                                                int categoryRatesIndex,
                                                final double[] inCategoryRates);

    public native int setTransitionMatrix(int instance, int matrixIndex, final double[] inMatrix, double paddedValue);

    public native int getTransitionMatrix(int instance, int matrixIndex, final double[] outMatrix);

	public native int convolveTransitionMatrices(int instance,
			                                     final int[] firstIndices, 
			                                     final int[] secondIndices,
			                                     final int[] resultIndices, 
			                                     int matrixCount);
    
    public native int updateTransitionMatrices(int instance, int eigenIndex,
                                               final int[] probabilityIndices,
                                               final int[] firstDerivativeIndices,
                                               final int[] secondDervativeIndices,
                                               final double[] edgeLengths,
                                               int count);

    public native int updateTransitionMatricesWithMultipleModels(
                                               int instance,
                                               final int[] eigenIndices,
                                               final int[] categoryRateIndices,
                                               final int[] probabilityIndices,
                                               final int[] firstDerivativeIndices,
                                               final int[] secondDervativeIndices,
                                               final double[] edgeLengths,
                                               int count);

    public native int updatePartials(final int instance,
                                     final int[] operations,
                                     int operationCount,
                                     int cumulativeScalingIndex);

    public native int updatePartialsByPartition(final int instance,
                                                final int[] operations,
                                                int operationCount);

    public native int waitForPartials(final int instance,
                                      final int[] destinationPartials,
                                      int destinationPartialsCount);

    public native int accumulateScaleFactors(final int instance,
                                             final int[] scaleIndices,
                                             final int count,
                                             final int cumulativeScalingIndex);

    public native int accumulateScaleFactorsByPartition(final int instance,
                                                        final int[] scaleIndices,
                                                        int count,
                                                        int cumulativeScaleIndex,
                                                        int partitionIndex);

    public native int removeScaleFactors(final int instance,
                                         final int[] scaleIndices,
                                         final int count,
                                         final int cumulativeScalingIndex);

    public native int removeScaleFactorsByPartition(final int instance,
                                         final int[] scaleIndices,
                                         final int count,
                                         final int cumulativeScalingIndex,
                                         final int partitionIndex);

    public native int resetScaleFactors(final int instance,
                                        final int cumulativeScalingIndex);

    public native int resetScaleFactorsByPartition(final int instance,
                                        final int cumulativeScalingIndex,
                                        final int partitionIndex);

    public native int copyScaleFactors(final int instance,
                                       final int destScalingIndex,
                                       final int srcScalingIndex);

    public native int calculateRootLogLikelihoods(int instance,
                                                  final int[] bufferIndices,
                                                  final int[] categoryWeightsIndices,
                                                  final int[] stateFrequenciesIndices,
                                                  final int[] cumulativeScaleIndices,
                                                  int count,
                                                  final double[] outSumLogLikelihood);

    public native int calculateRootLogLikelihoodsByPartition(int instance,
                                                  final int[] bufferIndices,
                                                  final int[] categoryWeightsIndices,
                                                  final int[] stateFrequenciesIndices,
                                                  final int[] cumulativeScaleIndices,
                                                  final int[] partitionIndices,
                                                  int partitionCount,
                                                  int count,
                                                  final double[] outSumLogLikelihoodByPartition,
                                                  final double[] outSumLogLikelihood);

    /*public native int calculateEdgeLogLikelihoods(int instance,
                                                  final int[] parentBufferIndices,
                                                  final int[] childBufferIndices,
                                                  final int[] probabilityIndices,
                                                  final int[] firstDerivativeIndices,
                                                  final int[] secondDerivativeIndices,
                                                  final int[] categoryWeightsIndices,
                                                  final int[] stateFrequenciesIndices,
                                                  final int[] scalingFactorsIndices,
                                                  int count,
                                                  final double[] outSumLogLikelihood,
                                                  final double[] outSumFirstDerivative,
                                                  final double[] outSumSecondDerivative);*/

    public native int getSiteLogLikelihoods(final int instance,
                                            final double[] outLogLikelihoods);

    /* Library loading routines */

    private static String getPlatformSpecificLibraryName()
    {
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        if (osName.startsWith("windows")) {
            if(osArch.equals("x86")||osArch.equals("i386")) return "hmsbeagle32";
            if(osArch.startsWith("amd64")||osArch.startsWith("x86_64")) return "hmsbeagle64";
        }
        return "hmsbeagle-jni";
    }

    public static void loadBeagleLibrary() throws UnsatisfiedLinkError {
        String path = "";
        if (System.getProperty("beagle.library.path") != null) {
            path = System.getProperty("beagle.library.path");
            if (path.length() > 0 && !path.endsWith("/")) {
                path += "/";
            }
        }

        System.loadLibrary(path + LIBRARY_NAME);
        INSTANCE = new BeagleJNIWrapper();
    }

    public static BeagleJNIWrapper INSTANCE;
}
