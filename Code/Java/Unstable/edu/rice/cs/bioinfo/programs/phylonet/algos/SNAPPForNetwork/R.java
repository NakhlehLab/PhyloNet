package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;


import org.apache.commons.math3.util.CombinatoricsUtils;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * This is a vector representing the number of each nucleotide base in the current branch.
 */
public class R {

    /**
     * The number of dims of the R vector.
     * Dims +1 = the number of nucleotides.
     * dims = 3 for the 4 nucleotide case.
     */
    private final static int dims = 3;


    /**
     * The total number of lineages.
     */
    int n;

    /**
     * The number of lineages with each nucleotide.
     */
    int[] values = new int[dims];



    int sum;


    /**
     * Creates an R vector with a certian number of lineages and lineage counts.
     * @param n The number of lineages.
     * @param values The lineage counts.
     */
    public R(int n, int[] values) {
        if (values.length != dims)
            throw new RuntimeException("Values is the wrong size.");

        if (!isValidValues(n,values))
            throw new RuntimeException("Invalid values.");

        this.n =n;
        this.values = values;

        this.sum = calculateSum(values);

    }


    private static int calculateSum(int[] arr)
    {
        int sum = 0;
        for (int elem : arr)
            sum += elem;

        return sum;
    }

    /**
     * Gets the size of a matrix used to store all R for 1<=n<=m.
     * @param m The maximum number of lineages.
     * @return The size of the matrix.
     */
    public static int getMatrixSize(int m){

        return calculateIndexedDimension(m,dims+1) -1;
    }

    /**
     * Gets the size of a matrix used to store all R with n total lineages.
     * @param n The total number of lineages.
     * @return The size of the matrix.
     */
    public static int getMatrixSizeWithoutN(int n) {

        return calculateIndexedDimension(n,dims);
    }



    @Override
    public String toString() {

        String sum = "R{";

        for (int type = 0; type < getNumberOfTypes(); type++)
        {
            sum += getTypeName(type) + "=" + getNum(type) + ",";
        }

        return  sum + '}';
    }


    /**
     * Gets the "type name" of a type. These names are arbitrary and are mainly just used for debugging.
     * @param type The type
     * @return The name of that type.
     */
    public String getTypeName(int type){
        switch (type)
        {
            case 0:
                return "A";
            case 1:
                return "C";
            case 2:
                return "T";
            case 3:
                return "G";

        }

        throw new RuntimeException("No such type");
    }


    /**
     * This calculates the number of indices for a vector of a given size and dimension.
     * How this works is not that clear.
     * @param index The size.
     * @param dimension The number of dimensions.
     * @return The number of indices.
     */
    private static int calculateIndexedDimension(int index, int dimension)
    {
        int result = 1;

        for (int i = 1;i <= dimension; i++)
        {
            result  *= (index + i);
        }

        for (int i = 1; i <= dimension; i++) {

            result /=i;
        }

        return result;

    }

    /**
     * Gets the index of this R in a matrix with only one n value.
     * @return The index.
     */
    public int getIndexWithoutN(){
        int total = 0;

        int totalCount = 0;
        for (int i = 0; i < dims; i++) {
            totalCount += values[i];
            total += calculateIndexedDimension(totalCount-1,i+1);
        }
        return total;
    }

    /**
     * Gets the index of this R in a matrix storing R's of varies n values.
     * @return The index.
     */
    public int getIndex(){
        if(n==0){
            return -1;
        }
        int total = calculateIndexedDimension(n-1,dims+1) -1;
        total += getIndexWithoutN();
        return total;
    }

    /**
     * Attempts to subtract one R from another.
     * Returns an optional as it can fail if it would result in negative size elements.
     * @param other The other R.
     * @return this - other
     */
    public R subtract(R other)
    {
        int newN = n - other.n;
        int[] newValues = subtractValues(values,other.values);
        if (isValidValues(newN,newValues))
            return new R(newN,newValues);
        else
            return null;
    }

    /**
     * Checks if the values are valid for a given n.
     * @param n The number of lineages.
     * @param values The lineage counts.
     * @return True if they are valid.
     */
    private static boolean isValidValues(int n, int[] values)
    {
        if (n < 0)
            return false;

        int remaining = n;

        for (int a : values)
        {
            if (a<0)
                return false;
            remaining-=a;
        }

        return remaining >=0;
    }

    /**
     * Calculates vector a-b.
     * @return A new array holding a-b.
     */
    private static int[] subtractValues(int[] a, int[]b)
    {
        int[] result = new int[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    /**
     * Returns the probability of selecting toSelect from this R.
     * Follows the multivariate hypergeometric distribution.
     * @param toSelect The R to select from.
     * @return The probablity of this selection.
     */
    public double getProbabilityOfSelecting(R toSelect)
    {
        double product = 1;
        for (int i =0; i < getNumberOfTypes(); i++)
        {
            product *= CombinatoricsUtils.binomialCoefficientDouble(getNum(i),toSelect.getNum(i));
        }
        product/= CombinatoricsUtils.binomialCoefficientDouble(n,toSelect.n);

        return product;
    }

    /**
     * Calculates a factor for selection probabilities.
     * @return n!/(v[0]! * v[1]! .. * v[numberOfTypes]!).
     */
    public double getLikelihoodProduct()
    {
        double numerator = CombinatoricsUtils.factorialDouble(n);

        for (int i = 0; i < getNumberOfTypes(); i++) {
            numerator/= CombinatoricsUtils.factorialDouble(getNum(i));
        }

        if (Double.isNaN(numerator)||Double.isInfinite(numerator))
            throw new RuntimeException("GetLikelihoodproduct is not returning a finite value.");

        return numerator;
    }


    /**
     * Gets the index for R in a square matrix of a size.
     * @param size The size of the square matrix.
     * @return The index.
     */
    public int getBoxIndex(int size)
    {
        int sum =n;

        for (int i = 0; i < dims; i++)
        {
            sum *= size;
            sum += getNum(i);
        }
        return sum;
    }

    /**
     * Returns an iterable over all values of R for a given n.
     * @param n The number of lineages.
     * @return An iterable to iterate over.
     */
    public static Iterable<R> loopOver(final int n) {
        return new Iterable<R>()
        {
            @Override
            public Iterator<R> iterator()
            {
                return new Iterator<R>()
                {

                    R last = null;

                    @Override
                    public void remove()
                    {
                        throw new UnsupportedOperationException();
                    }


                    @Override
                    public boolean hasNext()
                    {
                        return last == null || last.values[dims - 1] != n;
                    }

                    @Override
                    public R next()
                    {
                        if (last == null)
                            last = new R(n, new int[dims]);
                        else
                        {
                            tryAdvance(0);
                        }
                        last.sum = calculateSum(last.values);
                        return last;
                    }

                    private void tryAdvance(int i)
                    {
                        if (i == -1)
                            throw new RuntimeException("Iterator went over limit");


                        if (last.sum == n)
                        {
                            last.values[i] = 0;
                            last.sum = calculateSum(last.values);
                            tryAdvance(i + 1);
                        } else
                            last.values[i]++;

                    }
                };
            }
        };

    }

    /**
     * @return The number of types for this R.
     */
    public static int getNumberOfTypes()
    {
        return dims+1;
    }


    /**
     * Gets the number of lineages of a certain type from the R vector.
     * @param type The type
     * @return The number of lineages of that type.
     */
    public int getNum(int type)
    {
        if (type != dims)
            return values[type];

        else
            return n-sum;
    }


    /**
     * Assuming the given R value represents only one lineage, returns the type of that lineage.
     * @return The type
     */
    public int getType() {
        if (n!= 1)
            throw new RuntimeException("A type should only have one allele");


        for (int i = 0; i < getNumberOfTypes();i++)
        {
            if (getNum(i) == 1)
                return i;
        }

        throw new RuntimeException("Should be either red or green");
    }

    /**
     * Equivalent to transition(fromType,toType).getIndex();
     */
    public int transitionIndex(int fromType, int toType) {

        if (fromType != dims)
            values[fromType] -=1;

        if (toType != dims)
            values[toType] += 1;

        int result =  getIndex();

        if (fromType != dims)
            values[fromType] +=1;

        if (toType != dims)
            values[toType] -= 1;

        return result;
    }

    /**
     * Simulates a mutation from the fromType to the toType.
     * @param fromType The source lineage type.
     * @param toType The destination lineage type.
     * @return A new R representing after the mutation.
     */
    public R transition(int fromType, int toType) {

        int[] newValues = values.clone();

        if (fromType != dims)
            newValues[fromType] -=1;

        if (toType != dims)
            newValues[toType] += 1;

        return new R(n,newValues);
    }

    /**
     * Simulates coalescing two of the lineages into one.
     * @param coalesceType The type to coalesce.
     * @return A new R after coalesing.
     */
    public R coalesce(int coalesceType) {

        int[] newValues = values.clone();

        if (coalesceType != dims)
            newValues[coalesceType] -=1;

        return new R(n-1,newValues);
    }

    /**
     * Equivalent to split(splitType).getIndex()
     */
    public int splitIndex(int splitType) {

        if (splitType != dims)
            values[splitType] += 1;

        n+= 1;

        int result = getIndex();
        if (splitType != dims)
            values[splitType] -= 1;

        n-=1;

        return result;
    }

    /**
     * Simulates a lineage splitting into two. (IE, the opposite of coalescing).
     * @param splitType The type to split.
     * @return A new R after this split.
     */
    public R split(int splitType) {
        int[] newValues = values.clone();

        if (splitType != dims)
            newValues[splitType] += 1;

        return new R(n+1,newValues);
    }

    public static void main(String[] args)
    {

        Map<Integer,String> seenIndices = new HashMap<Integer,String>();

        for (int n = 1; n <=20 ;n++) {


            for (R r : R.loopOver(n)) {
                int firstIndex = r.getIndex();

                if (seenIndices.containsKey(firstIndex)) {
                    throw new RuntimeException("Fail");
                }

                seenIndices.put(firstIndex, r.toString());

                System.out.println(firstIndex + "," + r.getIndexWithoutN() + " , " + r);

            }
        }

    }

}
