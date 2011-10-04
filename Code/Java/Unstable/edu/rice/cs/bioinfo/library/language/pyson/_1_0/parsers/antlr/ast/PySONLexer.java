// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2011-10-03 18:32:48

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONLexer extends Lexer {
    public static final int EOF=-1;
    public static final int T__23=23;
    public static final int T__24=24;
    public static final int T__25=25;
    public static final int BEGIN=4;
    public static final int DEFAULT_INDICATOR=5;
    public static final int ELSE=6;
    public static final int END=7;
    public static final int ID=8;
    public static final int ID_SET=9;
    public static final int NESTED_ML_COMMENT=10;
    public static final int NETWORK=11;
    public static final int NETWORKS=12;
    public static final int PHYLONET=13;
    public static final int QUOTE=14;
    public static final int ROOTAGE_QUALIFIER=15;
    public static final int START=16;
    public static final int TAXON_SET_LIST=17;
    public static final int TRANSLATE=18;
    public static final int TREE=19;
    public static final int TREES=20;
    public static final int UTREE=21;
    public static final int WS=22;

    @Override   
    public void displayRecognitionError(String[] tokenNames, RecognitionException e) 
    {
    // dont show erors on std error
    }


    // delegates
    // delegators
    public Lexer[] getDelegates() {
        return new Lexer[] {};
    }

    public PySONLexer() {} 
    public PySONLexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public PySONLexer(CharStream input, RecognizerSharedState state) {
        super(input,state);
    }
    public String getGrammarFileName() { return "D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g"; }

    // $ANTLR start "T__23"
    public final void mT__23() throws RecognitionException {
        try {
            int _type = T__23;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:13:7: ( ',' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:13:9: ','
            {
            match(','); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__23"

    // $ANTLR start "T__24"
    public final void mT__24() throws RecognitionException {
        try {
            int _type = T__24;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:14:7: ( ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:14:9: ';'
            {
            match(';'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__24"

    // $ANTLR start "T__25"
    public final void mT__25() throws RecognitionException {
        try {
            int _type = T__25;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:15:7: ( '=' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:15:9: '='
            {
            match('='); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__25"

    // $ANTLR start "DEFAULT_INDICATOR"
    public final void mDEFAULT_INDICATOR() throws RecognitionException {
        try {
            int _type = DEFAULT_INDICATOR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:2: ( '*' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:4: '*'
            {
            match('*'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "DEFAULT_INDICATOR"

    // $ANTLR start "PHYLONET"
    public final void mPHYLONET() throws RecognitionException {
        try {
            int _type = PHYLONET;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:9: ( ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'Y' | 'y' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:11: ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'Y' | 'y' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' )
            {
            if ( input.LA(1)=='P'||input.LA(1)=='p' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='H'||input.LA(1)=='h' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='Y'||input.LA(1)=='y' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='L'||input.LA(1)=='l' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "PHYLONET"

    // $ANTLR start "TREE"
    public final void mTREE() throws RecognitionException {
        try {
            int _type = TREE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:124:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:124:10: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
            {
            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TREE"

    // $ANTLR start "UTREE"
    public final void mUTREE() throws RecognitionException {
        try {
            int _type = UTREE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:126:8: ( ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:126:11: ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
            {
            if ( input.LA(1)=='U'||input.LA(1)=='u' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "UTREE"

    // $ANTLR start "NETWORK"
    public final void mNETWORK() throws RecognitionException {
        try {
            int _type = NETWORK;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:128:9: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:128:11: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' )
            {
            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='W'||input.LA(1)=='w' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='K'||input.LA(1)=='k' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NETWORK"

    // $ANTLR start "START"
    public final void mSTART() throws RecognitionException {
        try {
            int _type = START;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:2: ( '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:4: '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' )
            {
            match('#'); 

            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='X'||input.LA(1)=='x' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='U'||input.LA(1)=='u' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "START"

    // $ANTLR start "BEGIN"
    public final void mBEGIN() throws RecognitionException {
        try {
            int _type = BEGIN;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:133:8: ( ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:133:10: ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' )
            {
            if ( input.LA(1)=='B'||input.LA(1)=='b' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='G'||input.LA(1)=='g' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='I'||input.LA(1)=='i' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "BEGIN"

    // $ANTLR start "TRANSLATE"
    public final void mTRANSLATE() throws RecognitionException {
        try {
            int _type = TRANSLATE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:136:2: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:136:4: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' )
            {
            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='L'||input.LA(1)=='l' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='A'||input.LA(1)=='a' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TRANSLATE"

    // $ANTLR start "NETWORKS"
    public final void mNETWORKS() throws RecognitionException {
        try {
            int _type = NETWORKS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:139:2: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:139:4: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' )
            {
            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='W'||input.LA(1)=='w' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='O'||input.LA(1)=='o' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='K'||input.LA(1)=='k' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NETWORKS"

    // $ANTLR start "TREES"
    public final void mTREES() throws RecognitionException {
        try {
            int _type = TREES;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:9: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' )
            {
            if ( input.LA(1)=='T'||input.LA(1)=='t' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='R'||input.LA(1)=='r' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='S'||input.LA(1)=='s' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TREES"

    // $ANTLR start "END"
    public final void mEND() throws RecognitionException {
        try {
            int _type = END;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:143:5: ( ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:143:7: ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' )
            {
            if ( input.LA(1)=='E'||input.LA(1)=='e' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='N'||input.LA(1)=='n' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            if ( input.LA(1)=='D'||input.LA(1)=='d' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "END"

    // $ANTLR start "ID"
    public final void mID() throws RecognitionException {
        try {
            int _type = ID;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:145:4: ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:145:6: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:145:6: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+
            int cnt1=0;
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( ((LA1_0 >= '0' && LA1_0 <= '9')||(LA1_0 >= 'A' && LA1_0 <= 'Z')||LA1_0=='_'||(LA1_0 >= 'a' && LA1_0 <= 'z')) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    if ( cnt1 >= 1 ) break loop1;
                        EarlyExitException eee =
                            new EarlyExitException(1, input);
                        throw eee;
                }
                cnt1++;
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ID"

    // $ANTLR start "QUOTE"
    public final void mQUOTE() throws RecognitionException {
        try {
            int _type = QUOTE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:147:7: ( '\"' (~ ( '\"' ) )* '\"' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:147:9: '\"' (~ ( '\"' ) )* '\"'
            {
            match('\"'); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:147:13: (~ ( '\"' ) )*
            loop2:
            do {
                int alt2=2;
                int LA2_0 = input.LA(1);

                if ( ((LA2_0 >= '\u0000' && LA2_0 <= '!')||(LA2_0 >= '#' && LA2_0 <= '\uFFFF')) ) {
                    alt2=1;
                }


                switch (alt2) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '!')||(input.LA(1) >= '#' && input.LA(1) <= '\uFFFF') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop2;
                }
            } while (true);


            match('\"'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "QUOTE"

    // $ANTLR start "ID_SET"
    public final void mID_SET() throws RecognitionException {
        try {
            int _type = ID_SET;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:2: ( '{' ( ( WS )* ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ ) ( WS )* ',' )* ( WS )* ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ ) ( WS )* '}' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:4: '{' ( ( WS )* ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ ) ( WS )* ',' )* ( WS )* ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ ) ( WS )* '}'
            {
            match('{'); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:8: ( ( WS )* ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ ) ( WS )* ',' )*
            loop6:
            do {
                int alt6=2;
                alt6 = dfa6.predict(input);
                switch (alt6) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:9: ( WS )* ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ ) ( WS )* ','
            	    {
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:9: ( WS )*
            	    loop3:
            	    do {
            	        int alt3=2;
            	        int LA3_0 = input.LA(1);

            	        if ( ((LA3_0 >= '\t' && LA3_0 <= '\n')||LA3_0=='\r'||LA3_0==' ') ) {
            	            alt3=1;
            	        }


            	        switch (alt3) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:9: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop3;
            	        }
            	    } while (true);


            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:13: ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ )
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:14: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+
            	    {
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:14: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+
            	    int cnt4=0;
            	    loop4:
            	    do {
            	        int alt4=2;
            	        int LA4_0 = input.LA(1);

            	        if ( ((LA4_0 >= '0' && LA4_0 <= '9')||(LA4_0 >= 'A' && LA4_0 <= 'Z')||LA4_0=='_'||(LA4_0 >= 'a' && LA4_0 <= 'z')) ) {
            	            alt4=1;
            	        }


            	        switch (alt4) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    	    {
            	    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
            	    	        input.consume();
            	    	    }
            	    	    else {
            	    	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	    	        recover(mse);
            	    	        throw mse;
            	    	    }


            	    	    }
            	    	    break;

            	    	default :
            	    	    if ( cnt4 >= 1 ) break loop4;
            	                EarlyExitException eee =
            	                    new EarlyExitException(4, input);
            	                throw eee;
            	        }
            	        cnt4++;
            	    } while (true);


            	    }


            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:55: ( WS )*
            	    loop5:
            	    do {
            	        int alt5=2;
            	        int LA5_0 = input.LA(1);

            	        if ( ((LA5_0 >= '\t' && LA5_0 <= '\n')||LA5_0=='\r'||LA5_0==' ') ) {
            	            alt5=1;
            	        }


            	        switch (alt5) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:55: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop5;
            	        }
            	    } while (true);


            	    match(','); 

            	    }
            	    break;

            	default :
            	    break loop6;
                }
            } while (true);


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:15: ( WS )*
            loop7:
            do {
                int alt7=2;
                int LA7_0 = input.LA(1);

                if ( ((LA7_0 >= '\t' && LA7_0 <= '\n')||LA7_0=='\r'||LA7_0==' ') ) {
                    alt7=1;
                }


                switch (alt7) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:15: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop7;
                }
            } while (true);


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:19: ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:20: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:20: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+
            int cnt8=0;
            loop8:
            do {
                int alt8=2;
                int LA8_0 = input.LA(1);

                if ( ((LA8_0 >= '0' && LA8_0 <= '9')||(LA8_0 >= 'A' && LA8_0 <= 'Z')||LA8_0=='_'||(LA8_0 >= 'a' && LA8_0 <= 'z')) ) {
                    alt8=1;
                }


                switch (alt8) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||input.LA(1)=='_'||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    if ( cnt8 >= 1 ) break loop8;
                        EarlyExitException eee =
                            new EarlyExitException(8, input);
                        throw eee;
                }
                cnt8++;
            } while (true);


            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:61: ( WS )*
            loop9:
            do {
                int alt9=2;
                int LA9_0 = input.LA(1);

                if ( ((LA9_0 >= '\t' && LA9_0 <= '\n')||LA9_0=='\r'||LA9_0==' ') ) {
                    alt9=1;
                }


                switch (alt9) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:61: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop9;
                }
            } while (true);


            match('}'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ID_SET"

    // $ANTLR start "TAXON_SET_LIST"
    public final void mTAXON_SET_LIST() throws RecognitionException {
        try {
            int _type = TAXON_SET_LIST;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:2: ( '(' ( ( WS )* ID_SET ( ( WS )* ',' )? )* ( ( WS )* ID_SET ( WS )* ) ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:4: '(' ( ( WS )* ID_SET ( ( WS )* ',' )? )* ( ( WS )* ID_SET ( WS )* ) ')'
            {
            match('('); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:8: ( ( WS )* ID_SET ( ( WS )* ',' )? )*
            loop13:
            do {
                int alt13=2;
                alt13 = dfa13.predict(input);
                switch (alt13) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:9: ( WS )* ID_SET ( ( WS )* ',' )?
            	    {
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:9: ( WS )*
            	    loop10:
            	    do {
            	        int alt10=2;
            	        int LA10_0 = input.LA(1);

            	        if ( ((LA10_0 >= '\t' && LA10_0 <= '\n')||LA10_0=='\r'||LA10_0==' ') ) {
            	            alt10=1;
            	        }


            	        switch (alt10) {
            	    	case 1 :
            	    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:9: WS
            	    	    {
            	    	    mWS(); 


            	    	    }
            	    	    break;

            	    	default :
            	    	    break loop10;
            	        }
            	    } while (true);


            	    mID_SET(); 


            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:20: ( ( WS )* ',' )?
            	    int alt12=2;
            	    alt12 = dfa12.predict(input);
            	    switch (alt12) {
            	        case 1 :
            	            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:21: ( WS )* ','
            	            {
            	            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:21: ( WS )*
            	            loop11:
            	            do {
            	                int alt11=2;
            	                int LA11_0 = input.LA(1);

            	                if ( ((LA11_0 >= '\t' && LA11_0 <= '\n')||LA11_0=='\r'||LA11_0==' ') ) {
            	                    alt11=1;
            	                }


            	                switch (alt11) {
            	            	case 1 :
            	            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:154:21: WS
            	            	    {
            	            	    mWS(); 


            	            	    }
            	            	    break;

            	            	default :
            	            	    break loop11;
            	                }
            	            } while (true);


            	            match(','); 

            	            }
            	            break;

            	    }


            	    }
            	    break;

            	default :
            	    break loop13;
                }
            } while (true);


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:14: ( ( WS )* ID_SET ( WS )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:15: ( WS )* ID_SET ( WS )*
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:15: ( WS )*
            loop14:
            do {
                int alt14=2;
                int LA14_0 = input.LA(1);

                if ( ((LA14_0 >= '\t' && LA14_0 <= '\n')||LA14_0=='\r'||LA14_0==' ') ) {
                    alt14=1;
                }


                switch (alt14) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:15: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop14;
                }
            } while (true);


            mID_SET(); 


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:26: ( WS )*
            loop15:
            do {
                int alt15=2;
                int LA15_0 = input.LA(1);

                if ( ((LA15_0 >= '\t' && LA15_0 <= '\n')||LA15_0=='\r'||LA15_0==' ') ) {
                    alt15=1;
                }


                switch (alt15) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:26: WS
            	    {
            	    mWS(); 


            	    }
            	    break;

            	default :
            	    break loop15;
                }
            } while (true);


            }


            match(')'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "TAXON_SET_LIST"

    // $ANTLR start "ROOTAGE_QUALIFIER"
    public final void mROOTAGE_QUALIFIER() throws RecognitionException {
        try {
            int _type = ROOTAGE_QUALIFIER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:158:2: ( '[' '&' ( 'r' | 'R' | 'u' | 'U' ) ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:158:4: '[' '&' ( 'r' | 'R' | 'u' | 'U' ) ']'
            {
            match('['); 

            match('&'); 

            if ( input.LA(1)=='R'||input.LA(1)=='U'||input.LA(1)=='r'||input.LA(1)=='u' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            match(']'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ROOTAGE_QUALIFIER"

    // $ANTLR start "NESTED_ML_COMMENT"
    public final void mNESTED_ML_COMMENT() throws RecognitionException {
        try {
            int _type = NESTED_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:3: ( '[' (~ ( '[' | ']' ) )* ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:6: '[' (~ ( '[' | ']' ) )* ']'
            {
            match('['); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:161:10: (~ ( '[' | ']' ) )*
            loop16:
            do {
                int alt16=2;
                int LA16_0 = input.LA(1);

                if ( ((LA16_0 >= '\u0000' && LA16_0 <= 'Z')||LA16_0=='\\'||(LA16_0 >= '^' && LA16_0 <= '\uFFFF')) ) {
                    alt16=1;
                }


                switch (alt16) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= 'Z')||input.LA(1)=='\\'||(input.LA(1) >= '^' && input.LA(1) <= '\uFFFF') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop16;
                }
            } while (true);


            match(']'); 

            _channel=HIDDEN;

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "NESTED_ML_COMMENT"

    // $ANTLR start "WS"
    public final void mWS() throws RecognitionException {
        try {
            int _type = WS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:164:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:164:9: ( ' ' | '\\t' | '\\r' | '\\n' )
            {
            if ( (input.LA(1) >= '\t' && input.LA(1) <= '\n')||input.LA(1)=='\r'||input.LA(1)==' ' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            _channel=HIDDEN;

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "WS"

    // $ANTLR start "ELSE"
    public final void mELSE() throws RecognitionException {
        try {
            int _type = ELSE;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:173:6: ( . )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:173:8: .
            {
            matchAny(); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "ELSE"

    public void mTokens() throws RecognitionException {
        // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:8: ( T__23 | T__24 | T__25 | DEFAULT_INDICATOR | PHYLONET | TREE | UTREE | NETWORK | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | ID | QUOTE | ID_SET | TAXON_SET_LIST | ROOTAGE_QUALIFIER | NESTED_ML_COMMENT | WS | ELSE )
        int alt17=22;
        alt17 = dfa17.predict(input);
        switch (alt17) {
            case 1 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:10: T__23
                {
                mT__23(); 


                }
                break;
            case 2 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:16: T__24
                {
                mT__24(); 


                }
                break;
            case 3 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:22: T__25
                {
                mT__25(); 


                }
                break;
            case 4 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:28: DEFAULT_INDICATOR
                {
                mDEFAULT_INDICATOR(); 


                }
                break;
            case 5 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:46: PHYLONET
                {
                mPHYLONET(); 


                }
                break;
            case 6 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:55: TREE
                {
                mTREE(); 


                }
                break;
            case 7 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:60: UTREE
                {
                mUTREE(); 


                }
                break;
            case 8 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:66: NETWORK
                {
                mNETWORK(); 


                }
                break;
            case 9 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:74: START
                {
                mSTART(); 


                }
                break;
            case 10 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:80: BEGIN
                {
                mBEGIN(); 


                }
                break;
            case 11 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:86: TRANSLATE
                {
                mTRANSLATE(); 


                }
                break;
            case 12 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:96: NETWORKS
                {
                mNETWORKS(); 


                }
                break;
            case 13 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:105: TREES
                {
                mTREES(); 


                }
                break;
            case 14 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:111: END
                {
                mEND(); 


                }
                break;
            case 15 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:115: ID
                {
                mID(); 


                }
                break;
            case 16 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:118: QUOTE
                {
                mQUOTE(); 


                }
                break;
            case 17 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:124: ID_SET
                {
                mID_SET(); 


                }
                break;
            case 18 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:131: TAXON_SET_LIST
                {
                mTAXON_SET_LIST(); 


                }
                break;
            case 19 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:146: ROOTAGE_QUALIFIER
                {
                mROOTAGE_QUALIFIER(); 


                }
                break;
            case 20 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:164: NESTED_ML_COMMENT
                {
                mNESTED_ML_COMMENT(); 


                }
                break;
            case 21 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:182: WS
                {
                mWS(); 


                }
                break;
            case 22 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:185: ELSE
                {
                mELSE(); 


                }
                break;

        }

    }


    protected DFA6 dfa6 = new DFA6(this);
    protected DFA13 dfa13 = new DFA13(this);
    protected DFA12 dfa12 = new DFA12(this);
    protected DFA17 dfa17 = new DFA17(this);
    static final String DFA6_eotS =
        "\6\uffff";
    static final String DFA6_eofS =
        "\6\uffff";
    static final String DFA6_minS =
        "\4\11\2\uffff";
    static final String DFA6_maxS =
        "\2\172\2\175\2\uffff";
    static final String DFA6_acceptS =
        "\4\uffff\1\2\1\1";
    static final String DFA6_specialS =
        "\6\uffff}>";
    static final String[] DFA6_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\17\uffff\12\2\7\uffff\32\2\4\uffff"+
            "\1\2\1\uffff\32\2",
            "\2\1\2\uffff\1\1\22\uffff\1\1\17\uffff\12\2\7\uffff\32\2\4"+
            "\uffff\1\2\1\uffff\32\2",
            "\2\3\2\uffff\1\3\22\uffff\1\3\13\uffff\1\5\3\uffff\12\2\7\uffff"+
            "\32\2\4\uffff\1\2\1\uffff\32\2\2\uffff\1\4",
            "\2\3\2\uffff\1\3\22\uffff\1\3\13\uffff\1\5\120\uffff\1\4",
            "",
            ""
    };

    static final short[] DFA6_eot = DFA.unpackEncodedString(DFA6_eotS);
    static final short[] DFA6_eof = DFA.unpackEncodedString(DFA6_eofS);
    static final char[] DFA6_min = DFA.unpackEncodedStringToUnsignedChars(DFA6_minS);
    static final char[] DFA6_max = DFA.unpackEncodedStringToUnsignedChars(DFA6_maxS);
    static final short[] DFA6_accept = DFA.unpackEncodedString(DFA6_acceptS);
    static final short[] DFA6_special = DFA.unpackEncodedString(DFA6_specialS);
    static final short[][] DFA6_transition;

    static {
        int numStates = DFA6_transitionS.length;
        DFA6_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA6_transition[i] = DFA.unpackEncodedString(DFA6_transitionS[i]);
        }
    }

    class DFA6 extends DFA {

        public DFA6(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 6;
            this.eot = DFA6_eot;
            this.eof = DFA6_eof;
            this.min = DFA6_min;
            this.max = DFA6_max;
            this.accept = DFA6_accept;
            this.special = DFA6_special;
            this.transition = DFA6_transition;
        }
        public String getDescription() {
            return "()* loopback of 150:8: ( ( WS )* ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) | '_' )+ ) ( WS )* ',' )*";
        }
    }
    static final String DFA13_eotS =
        "\13\uffff";
    static final String DFA13_eofS =
        "\13\uffff";
    static final String DFA13_minS =
        "\11\11\2\uffff";
    static final String DFA13_maxS =
        "\2\173\2\172\2\175\1\172\2\173\2\uffff";
    static final String DFA13_acceptS =
        "\11\uffff\1\2\1\1";
    static final String DFA13_specialS =
        "\13\uffff}>";
    static final String[] DFA13_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\132\uffff\1\2",
            "\2\1\2\uffff\1\1\22\uffff\1\1\132\uffff\1\2",
            "\2\3\2\uffff\1\3\22\uffff\1\3\17\uffff\12\4\7\uffff\32\4\4"+
            "\uffff\1\4\1\uffff\32\4",
            "\2\3\2\uffff\1\3\22\uffff\1\3\17\uffff\12\4\7\uffff\32\4\4"+
            "\uffff\1\4\1\uffff\32\4",
            "\2\5\2\uffff\1\5\22\uffff\1\5\13\uffff\1\6\3\uffff\12\4\7\uffff"+
            "\32\4\4\uffff\1\4\1\uffff\32\4\2\uffff\1\7",
            "\2\5\2\uffff\1\5\22\uffff\1\5\13\uffff\1\6\120\uffff\1\7",
            "\2\3\2\uffff\1\3\22\uffff\1\3\17\uffff\12\4\7\uffff\32\4\4"+
            "\uffff\1\4\1\uffff\32\4",
            "\2\10\2\uffff\1\10\22\uffff\1\10\10\uffff\1\11\2\uffff\1\12"+
            "\116\uffff\1\12",
            "\2\10\2\uffff\1\10\22\uffff\1\10\10\uffff\1\11\2\uffff\1\12"+
            "\116\uffff\1\12",
            "",
            ""
    };

    static final short[] DFA13_eot = DFA.unpackEncodedString(DFA13_eotS);
    static final short[] DFA13_eof = DFA.unpackEncodedString(DFA13_eofS);
    static final char[] DFA13_min = DFA.unpackEncodedStringToUnsignedChars(DFA13_minS);
    static final char[] DFA13_max = DFA.unpackEncodedStringToUnsignedChars(DFA13_maxS);
    static final short[] DFA13_accept = DFA.unpackEncodedString(DFA13_acceptS);
    static final short[] DFA13_special = DFA.unpackEncodedString(DFA13_specialS);
    static final short[][] DFA13_transition;

    static {
        int numStates = DFA13_transitionS.length;
        DFA13_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA13_transition[i] = DFA.unpackEncodedString(DFA13_transitionS[i]);
        }
    }

    class DFA13 extends DFA {

        public DFA13(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 13;
            this.eot = DFA13_eot;
            this.eof = DFA13_eof;
            this.min = DFA13_min;
            this.max = DFA13_max;
            this.accept = DFA13_accept;
            this.special = DFA13_special;
            this.transition = DFA13_transition;
        }
        public String getDescription() {
            return "()* loopback of 154:8: ( ( WS )* ID_SET ( ( WS )* ',' )? )*";
        }
    }
    static final String DFA12_eotS =
        "\4\uffff";
    static final String DFA12_eofS =
        "\4\uffff";
    static final String DFA12_minS =
        "\2\11\2\uffff";
    static final String DFA12_maxS =
        "\2\173\2\uffff";
    static final String DFA12_acceptS =
        "\2\uffff\1\1\1\2";
    static final String DFA12_specialS =
        "\4\uffff}>";
    static final String[] DFA12_transitionS = {
            "\2\1\2\uffff\1\1\22\uffff\1\1\13\uffff\1\2\116\uffff\1\3",
            "\2\1\2\uffff\1\1\22\uffff\1\1\13\uffff\1\2\116\uffff\1\3",
            "",
            ""
    };

    static final short[] DFA12_eot = DFA.unpackEncodedString(DFA12_eotS);
    static final short[] DFA12_eof = DFA.unpackEncodedString(DFA12_eofS);
    static final char[] DFA12_min = DFA.unpackEncodedStringToUnsignedChars(DFA12_minS);
    static final char[] DFA12_max = DFA.unpackEncodedStringToUnsignedChars(DFA12_maxS);
    static final short[] DFA12_accept = DFA.unpackEncodedString(DFA12_acceptS);
    static final short[] DFA12_special = DFA.unpackEncodedString(DFA12_specialS);
    static final short[][] DFA12_transition;

    static {
        int numStates = DFA12_transitionS.length;
        DFA12_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA12_transition[i] = DFA.unpackEncodedString(DFA12_transitionS[i]);
        }
    }

    class DFA12 extends DFA {

        public DFA12(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 12;
            this.eot = DFA12_eot;
            this.eof = DFA12_eof;
            this.min = DFA12_min;
            this.max = DFA12_max;
            this.accept = DFA12_accept;
            this.special = DFA12_special;
            this.transition = DFA12_transition;
        }
        public String getDescription() {
            return "154:20: ( ( WS )* ',' )?";
        }
    }
    static final String DFA17_eotS =
        "\5\uffff\4\30\1\22\2\30\1\uffff\4\22\6\uffff\1\30\1\uffff\3\30\1"+
        "\uffff\2\30\6\uffff\6\30\1\63\1\uffff\1\30\1\67\4\30\2\uffff\1\30"+
        "\1\76\1\uffff\1\30\1\100\1\30\1\102\1\uffff\1\30\1\uffff\1\30\1"+
        "\uffff\1\30\1\uffff\2\30\1\111\1\112\1\30\1\114\2\uffff\1\115\2"+
        "\uffff";
    static final String DFA17_eofS =
        "\116\uffff";
    static final String DFA17_minS =
        "\1\0\4\uffff\1\110\1\122\1\124\1\105\1\116\1\105\1\116\1\uffff\1"+
        "\0\2\11\1\0\6\uffff\1\131\1\uffff\1\101\1\122\1\124\1\uffff\1\107"+
        "\1\104\3\uffff\1\0\2\uffff\1\114\1\105\1\116\1\105\1\127\1\111\1"+
        "\60\1\0\1\117\1\60\1\123\1\105\1\117\1\116\2\uffff\1\116\1\60\1"+
        "\uffff\1\114\1\60\1\122\1\60\1\uffff\1\105\1\uffff\1\101\1\uffff"+
        "\1\113\1\uffff\2\124\2\60\1\105\1\60\2\uffff\1\60\2\uffff";
    static final String DFA17_maxS =
        "\1\uffff\4\uffff\1\150\1\162\1\164\1\145\1\156\1\145\1\156\1\uffff"+
        "\1\uffff\1\172\1\173\1\uffff\6\uffff\1\171\1\uffff\1\145\1\162\1"+
        "\164\1\uffff\1\147\1\144\3\uffff\1\uffff\2\uffff\1\154\1\145\1\156"+
        "\1\145\1\167\1\151\1\172\1\uffff\1\157\1\172\1\163\1\145\1\157\1"+
        "\156\2\uffff\1\156\1\172\1\uffff\1\154\1\172\1\162\1\172\1\uffff"+
        "\1\145\1\uffff\1\141\1\uffff\1\153\1\uffff\2\164\2\172\1\145\1\172"+
        "\2\uffff\1\172\2\uffff";
    static final String DFA17_acceptS =
        "\1\uffff\1\1\1\2\1\3\1\4\7\uffff\1\17\4\uffff\1\25\1\26\1\1\1\2"+
        "\1\3\1\4\1\uffff\1\17\3\uffff\1\11\2\uffff\1\20\1\21\1\22\1\uffff"+
        "\1\24\1\25\16\uffff\1\16\1\23\2\uffff\1\6\4\uffff\1\23\1\uffff\1"+
        "\15\1\uffff\1\7\1\uffff\1\12\6\uffff\1\10\1\5\1\uffff\1\14\1\13";
    static final String DFA17_specialS =
        "\1\4\14\uffff\1\3\2\uffff\1\2\21\uffff\1\0\11\uffff\1\1\41\uffff}>";
    static final String[] DFA17_transitionS = {
            "\11\22\2\21\2\22\1\21\22\22\1\21\1\22\1\15\1\11\4\22\1\17\1"+
            "\22\1\4\1\22\1\1\3\22\12\14\1\22\1\2\1\22\1\3\3\22\1\14\1\12"+
            "\2\14\1\13\10\14\1\10\1\14\1\5\3\14\1\6\1\7\5\14\1\20\3\22\1"+
            "\14\1\22\1\14\1\12\2\14\1\13\10\14\1\10\1\14\1\5\3\14\1\6\1"+
            "\7\5\14\1\16\uff84\22",
            "",
            "",
            "",
            "",
            "\1\27\37\uffff\1\27",
            "\1\31\37\uffff\1\31",
            "\1\32\37\uffff\1\32",
            "\1\33\37\uffff\1\33",
            "\1\34\37\uffff\1\34",
            "\1\35\37\uffff\1\35",
            "\1\36\37\uffff\1\36",
            "",
            "\0\37",
            "\2\40\2\uffff\1\40\22\uffff\1\40\17\uffff\12\40\7\uffff\32"+
            "\40\4\uffff\1\40\1\uffff\32\40",
            "\2\41\2\uffff\1\41\22\uffff\1\41\132\uffff\1\41",
            "\46\43\1\42\64\43\1\uffff\uffa4\43",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\45\37\uffff\1\45",
            "",
            "\1\47\3\uffff\1\46\33\uffff\1\47\3\uffff\1\46",
            "\1\50\37\uffff\1\50",
            "\1\51\37\uffff\1\51",
            "",
            "\1\52\37\uffff\1\52",
            "\1\53\37\uffff\1\53",
            "",
            "",
            "",
            "\122\43\1\54\2\43\1\54\5\43\1\uffff\26\43\1\54\2\43\1\54\uff8a"+
            "\43",
            "",
            "",
            "\1\55\37\uffff\1\55",
            "\1\56\37\uffff\1\56",
            "\1\57\37\uffff\1\57",
            "\1\60\37\uffff\1\60",
            "\1\61\37\uffff\1\61",
            "\1\62\37\uffff\1\62",
            "\12\30\7\uffff\32\30\4\uffff\1\30\1\uffff\32\30",
            "\133\43\1\uffff\1\43\1\64\uffa2\43",
            "\1\65\37\uffff\1\65",
            "\12\30\7\uffff\22\30\1\66\7\30\4\uffff\1\30\1\uffff\22\30\1"+
            "\66\7\30",
            "\1\70\37\uffff\1\70",
            "\1\71\37\uffff\1\71",
            "\1\72\37\uffff\1\72",
            "\1\73\37\uffff\1\73",
            "",
            "",
            "\1\75\37\uffff\1\75",
            "\12\30\7\uffff\32\30\4\uffff\1\30\1\uffff\32\30",
            "",
            "\1\77\37\uffff\1\77",
            "\12\30\7\uffff\32\30\4\uffff\1\30\1\uffff\32\30",
            "\1\101\37\uffff\1\101",
            "\12\30\7\uffff\32\30\4\uffff\1\30\1\uffff\32\30",
            "",
            "\1\103\37\uffff\1\103",
            "",
            "\1\104\37\uffff\1\104",
            "",
            "\1\105\37\uffff\1\105",
            "",
            "\1\106\37\uffff\1\106",
            "\1\107\37\uffff\1\107",
            "\12\30\7\uffff\22\30\1\110\7\30\4\uffff\1\30\1\uffff\22\30"+
            "\1\110\7\30",
            "\12\30\7\uffff\32\30\4\uffff\1\30\1\uffff\32\30",
            "\1\113\37\uffff\1\113",
            "\12\30\7\uffff\32\30\4\uffff\1\30\1\uffff\32\30",
            "",
            "",
            "\12\30\7\uffff\32\30\4\uffff\1\30\1\uffff\32\30",
            "",
            ""
    };

    static final short[] DFA17_eot = DFA.unpackEncodedString(DFA17_eotS);
    static final short[] DFA17_eof = DFA.unpackEncodedString(DFA17_eofS);
    static final char[] DFA17_min = DFA.unpackEncodedStringToUnsignedChars(DFA17_minS);
    static final char[] DFA17_max = DFA.unpackEncodedStringToUnsignedChars(DFA17_maxS);
    static final short[] DFA17_accept = DFA.unpackEncodedString(DFA17_acceptS);
    static final short[] DFA17_special = DFA.unpackEncodedString(DFA17_specialS);
    static final short[][] DFA17_transition;

    static {
        int numStates = DFA17_transitionS.length;
        DFA17_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA17_transition[i] = DFA.unpackEncodedString(DFA17_transitionS[i]);
        }
    }

    class DFA17 extends DFA {

        public DFA17(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 17;
            this.eot = DFA17_eot;
            this.eof = DFA17_eof;
            this.min = DFA17_min;
            this.max = DFA17_max;
            this.accept = DFA17_accept;
            this.special = DFA17_special;
            this.transition = DFA17_transition;
        }
        public String getDescription() {
            return "1:1: Tokens : ( T__23 | T__24 | T__25 | DEFAULT_INDICATOR | PHYLONET | TREE | UTREE | NETWORK | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | ID | QUOTE | ID_SET | TAXON_SET_LIST | ROOTAGE_QUALIFIER | NESTED_ML_COMMENT | WS | ELSE );";
        }
        public int specialStateTransition(int s, IntStream _input) throws NoViableAltException {
            IntStream input = _input;
        	int _s = s;
            switch ( s ) {
                    case 0 : 
                        int LA17_34 = input.LA(1);

                        s = -1;
                        if ( (LA17_34=='R'||LA17_34=='U'||LA17_34=='r'||LA17_34=='u') ) {s = 44;}

                        else if ( ((LA17_34 >= '\u0000' && LA17_34 <= 'Q')||(LA17_34 >= 'S' && LA17_34 <= 'T')||(LA17_34 >= 'V' && LA17_34 <= 'Z')||(LA17_34 >= '\\' && LA17_34 <= 'q')||(LA17_34 >= 's' && LA17_34 <= 't')||(LA17_34 >= 'v' && LA17_34 <= '\uFFFF')) ) {s = 35;}

                        if ( s>=0 ) return s;
                        break;

                    case 1 : 
                        int LA17_44 = input.LA(1);

                        s = -1;
                        if ( (LA17_44==']') ) {s = 52;}

                        else if ( ((LA17_44 >= '\u0000' && LA17_44 <= 'Z')||LA17_44=='\\'||(LA17_44 >= '^' && LA17_44 <= '\uFFFF')) ) {s = 35;}

                        if ( s>=0 ) return s;
                        break;

                    case 2 : 
                        int LA17_16 = input.LA(1);

                        s = -1;
                        if ( (LA17_16=='&') ) {s = 34;}

                        else if ( ((LA17_16 >= '\u0000' && LA17_16 <= '%')||(LA17_16 >= '\'' && LA17_16 <= 'Z')||(LA17_16 >= '\\' && LA17_16 <= '\uFFFF')) ) {s = 35;}

                        else s = 18;

                        if ( s>=0 ) return s;
                        break;

                    case 3 : 
                        int LA17_13 = input.LA(1);

                        s = -1;
                        if ( ((LA17_13 >= '\u0000' && LA17_13 <= '\uFFFF')) ) {s = 31;}

                        else s = 18;

                        if ( s>=0 ) return s;
                        break;

                    case 4 : 
                        int LA17_0 = input.LA(1);

                        s = -1;
                        if ( (LA17_0==',') ) {s = 1;}

                        else if ( (LA17_0==';') ) {s = 2;}

                        else if ( (LA17_0=='=') ) {s = 3;}

                        else if ( (LA17_0=='*') ) {s = 4;}

                        else if ( (LA17_0=='P'||LA17_0=='p') ) {s = 5;}

                        else if ( (LA17_0=='T'||LA17_0=='t') ) {s = 6;}

                        else if ( (LA17_0=='U'||LA17_0=='u') ) {s = 7;}

                        else if ( (LA17_0=='N'||LA17_0=='n') ) {s = 8;}

                        else if ( (LA17_0=='#') ) {s = 9;}

                        else if ( (LA17_0=='B'||LA17_0=='b') ) {s = 10;}

                        else if ( (LA17_0=='E'||LA17_0=='e') ) {s = 11;}

                        else if ( ((LA17_0 >= '0' && LA17_0 <= '9')||LA17_0=='A'||(LA17_0 >= 'C' && LA17_0 <= 'D')||(LA17_0 >= 'F' && LA17_0 <= 'M')||LA17_0=='O'||(LA17_0 >= 'Q' && LA17_0 <= 'S')||(LA17_0 >= 'V' && LA17_0 <= 'Z')||LA17_0=='_'||LA17_0=='a'||(LA17_0 >= 'c' && LA17_0 <= 'd')||(LA17_0 >= 'f' && LA17_0 <= 'm')||LA17_0=='o'||(LA17_0 >= 'q' && LA17_0 <= 's')||(LA17_0 >= 'v' && LA17_0 <= 'z')) ) {s = 12;}

                        else if ( (LA17_0=='\"') ) {s = 13;}

                        else if ( (LA17_0=='{') ) {s = 14;}

                        else if ( (LA17_0=='(') ) {s = 15;}

                        else if ( (LA17_0=='[') ) {s = 16;}

                        else if ( ((LA17_0 >= '\t' && LA17_0 <= '\n')||LA17_0=='\r'||LA17_0==' ') ) {s = 17;}

                        else if ( ((LA17_0 >= '\u0000' && LA17_0 <= '\b')||(LA17_0 >= '\u000B' && LA17_0 <= '\f')||(LA17_0 >= '\u000E' && LA17_0 <= '\u001F')||LA17_0=='!'||(LA17_0 >= '$' && LA17_0 <= '\'')||LA17_0==')'||LA17_0=='+'||(LA17_0 >= '-' && LA17_0 <= '/')||LA17_0==':'||LA17_0=='<'||(LA17_0 >= '>' && LA17_0 <= '@')||(LA17_0 >= '\\' && LA17_0 <= '^')||LA17_0=='`'||(LA17_0 >= '|' && LA17_0 <= '\uFFFF')) ) {s = 18;}

                        if ( s>=0 ) return s;
                        break;
            }
            NoViableAltException nvae =
                new NoViableAltException(getDescription(), 17, _s, input);
            error(nvae);
            throw nvae;
        }

    }
 

}