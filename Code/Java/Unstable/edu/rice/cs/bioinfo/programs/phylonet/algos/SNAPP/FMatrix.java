package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP;

import jeigen.DenseMatrix;

public class FMatrix
{
    double [] arr;
    public FMatrix(int mx)
    {
        arr = new double[R.getMatrixSize(mx)];

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
        return arr[r.getIndex()];
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

    public String toString(){
        String exp = "";
        for(int i=0; i<arr.length; i++){
            exp += arr[i] + "\n";
        }
        return exp;
    }
}
