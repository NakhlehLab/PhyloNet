package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import jeigen.DenseMatrix;

public class FMatrix
{
    double [] arr;
    int mx;
    boolean hasEmptyR;

    public FMatrix(int mxValue, boolean emptyRValue)
    {
        mx = mxValue;
        arr = new double[R.getMatrixSize(mx)];
        hasEmptyR = emptyRValue;
    }



    public FMatrix(int mxValue, double[] arrValue, boolean emptyRValue)
    {
        mx = mxValue;
        arr = arrValue;
        hasEmptyR = emptyRValue;
    }


    public boolean isArrAllZero(){
        for(double value: arr){
            if(value != 0){
                return false;
            }
        }
        return true;
    }

    public DenseMatrix getMatrix(){
        DenseMatrix result = new DenseMatrix(arr.length,1);
        for (int i = 0; i < arr.length; i++) {
            result.set(i,0,arr[i]);
        }
        return result;
    }

    public double get(R r)
    {
        int index = r.getIndex();
        if(index == -1){
            if(hasEmptyR){
                return 1;
            }
            else{
                return 0;
            }
        }
        else{
            return arr[r.getIndex()];
        }
    }

    public void setEmptyR(boolean has){
        hasEmptyR = has;
    }

    public boolean ifHasEmptyR(){
        return hasEmptyR;
    }


    public void set(R r, double val)
    {
        arr[r.getIndex()] = val;
    }

    public void setMatrix(DenseMatrix matrix) {

        for (int i = 0; i < matrix.rows; i++)
        {
            arr[i] = matrix.get(i,0);
        }
    }

    public double[] getArr(){
        return arr;
    }

    /*
    public String toString(){
        String exp = "";
        for(int i=0; i<arr.length; i++){
            exp += arr[i] + "\n";
        }
        return exp;
    }
    */

    public String toString(){
        String exp = "";
        for (int n = 1; n <= mx; n++) {
            for (R r : R.loopOver(n)) {
                exp += r + ": " + arr[r.getIndex()] + "\n";
            }
        }
        return exp;
    }
}
