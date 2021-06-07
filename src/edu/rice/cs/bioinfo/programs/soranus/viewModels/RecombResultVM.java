package edu.rice.cs.bioinfo.programs.soranus.viewModels;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/25/13
 * Time: 4:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class RecombResultVM implements DocumentVM
{
    public final String Sequence1Label;

    public final String Sequence2Label;

    public final String Sequence3Label;

    public final String Sequence4Label;

    public final String Site1Label;

    public final String Site2Label;

    public final char Sequence0Site0;

    public final char Sequence1Site0;

    public final char Sequence2Site0;

    public final char Sequence3Site0;

    public final char Sequence0Site1;

    public final char Sequence1Site1;

    public final char Sequence2Site1;

    public final char Sequence3Site1;

    public RecombResultVM(String sequence1Label, String sequence2Label, String sequence3Label, String sequence4Label, String site1Label, String site2Label, char sequence0Site0, char sequence1Site0, char sequence2Site0, char sequence3Site0, char sequence0Site1,
                          char sequence1Site1, char sequence2Site1, char sequence3Site1)
    {
        Sequence1Label = sequence1Label;
        Sequence2Label = sequence2Label;
        Sequence3Label = sequence3Label;
        Sequence4Label = sequence4Label;
        Site1Label = site1Label;
        Site2Label = site2Label;
        Sequence0Site0 = sequence0Site0;
        Sequence1Site0 = sequence1Site0;
        Sequence2Site0 = sequence2Site0;
        Sequence3Site0 = sequence3Site0;
        Sequence0Site1 = sequence0Site1;
        Sequence1Site1 = sequence1Site1;
        Sequence2Site1 = sequence2Site1;
        Sequence3Site1 = sequence3Site1;
    }

    public <R, E extends Exception> R execute(DocumentVMAlgo<R, E> algo) throws E
    {
        return algo.forRecombResult(this);
    }
}
