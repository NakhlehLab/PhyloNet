package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.datatype;

import java.util.Arrays;
import java.util.List;

/**
 * Created by wendingqiao on 5/3/16.
 */
public class Nucleotide extends DataType.Base {

    /*int[][] x = {
            {0},  // A
            {1},  // C
            {2},  // G
            {3},  // T
            {3},  // U
            {0, 2}, // R
            {1, 3}, // Y
            {0, 1}, // M
            {0, 3}, // W
            {1, 2}, // S
            {2, 3}, // K
            {1, 2, 3}, // B
            {0, 2, 3}, // D
            {0, 1, 3}, // H
            {0, 1, 2}, // V
            {0, 1, 2, 3}, // N
            {0, 1, 2, 3}, // X
            {0, 1, 2, 3}, // -
            {0, 1, 2, 3}, // ?
    };*/

    int[][] x = {
            {0},  // 0
            {1}   // 1
    };

    public Nucleotide() {
        /*stateCount = 4;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "ACGTURYMWSKBDHVNX" + GAP_CHAR + MISSING_CHAR;*/
        stateCount = 2;
        mapCodeToStateSet = x;
        codeLength = 1;
        //codeMap = "01";
        codeMap = "0123456789" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "nucleotide";
    }

    // test
    public static void main(String[] args) {
        Nucleotide nc = new Nucleotide();
        List<Integer> list = nc.stringToState("ACGT-");
        boolean correct = list.get(4) == 17;
        System.out.println(correct);
        for(int i = 0; i < 4; i++) {
            correct |= list.get(i) == i;
        }
        System.out.println(correct);
        System.out.println(Arrays.toString(nc.getStateSet(0))); // TFFF
        System.out.println(Arrays.toString(nc.getStatesForCode(17))); // 0123
        System.out.println((nc.getChar(0) == 'A') && (nc.getChar(2) == 'C'));
    }

}