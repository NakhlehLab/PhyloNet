// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2011-09-07 17:34:43

package edu.rice.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONLexer extends Lexer {
    public static final int EOF=-1;
    public static final int T__20=20;
    public static final int T__21=21;
    public static final int T__22=22;
    public static final int BEGIN=4;
    public static final int DEFAULT_INDICATOR=5;
    public static final int ELSE=6;
    public static final int END=7;
    public static final int ID=8;
    public static final int NESTED_ML_COMMENT=9;
    public static final int NETWORK=10;
    public static final int NETWORKS=11;
    public static final int PHYLONET=12;
    public static final int QUOTE=13;
    public static final int START=14;
    public static final int TRANSLATE=15;
    public static final int TREE=16;
    public static final int TREES=17;
    public static final int UTREE=18;
    public static final int WS=19;

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

    // $ANTLR start "T__20"
    public final void mT__20() throws RecognitionException {
        try {
            int _type = T__20;
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
    // $ANTLR end "T__20"

    // $ANTLR start "T__21"
    public final void mT__21() throws RecognitionException {
        try {
            int _type = T__21;
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
    // $ANTLR end "T__21"

    // $ANTLR start "T__22"
    public final void mT__22() throws RecognitionException {
        try {
            int _type = T__22;
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
    // $ANTLR end "T__22"

    // $ANTLR start "DEFAULT_INDICATOR"
    public final void mDEFAULT_INDICATOR() throws RecognitionException {
        try {
            int _type = DEFAULT_INDICATOR;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:116:2: ( '*' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:116:4: '*'
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:118:9: ( ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'Y' | 'y' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:118:11: ( 'P' | 'p' ) ( 'H' | 'h' ) ( 'Y' | 'y' ) ( 'L' | 'l' ) ( 'O' | 'o' ) ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:10: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:8: ( ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:11: ( 'U' | 'u' ) ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:124:9: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:124:11: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:2: ( '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:4: '#' ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'X' | 'x' ) ( 'U' | 'u' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:129:8: ( ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:129:10: ( 'B' | 'b' ) ( 'E' | 'e' ) ( 'G' | 'g' ) ( 'I' | 'i' ) ( 'N' | 'n' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:132:2: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:132:4: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'A' | 'a' ) ( 'N' | 'n' ) ( 'S' | 's' ) ( 'L' | 'l' ) ( 'A' | 'a' ) ( 'T' | 't' ) ( 'E' | 'e' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:135:2: ( ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:135:4: ( 'N' | 'n' ) ( 'E' | 'e' ) ( 'T' | 't' ) ( 'W' | 'w' ) ( 'O' | 'o' ) ( 'R' | 'r' ) ( 'K' | 'k' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:137:7: ( ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:137:9: ( 'T' | 't' ) ( 'R' | 'r' ) ( 'E' | 'e' ) ( 'E' | 'e' ) ( 'S' | 's' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:139:5: ( ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:139:7: ( 'E' | 'e' ) ( 'N' | 'n' ) ( 'D' | 'd' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:4: ( ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) )+ )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:6: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) )+
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:6: ( ( 'a' .. 'z' ) | ( 'A' .. 'Z' ) | ( '0' .. '9' ) )+
            int cnt1=0;
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( ((LA1_0 >= '0' && LA1_0 <= '9')||(LA1_0 >= 'A' && LA1_0 <= 'Z')||(LA1_0 >= 'a' && LA1_0 <= 'z')) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9')||(input.LA(1) >= 'A' && input.LA(1) <= 'Z')||(input.LA(1) >= 'a' && input.LA(1) <= 'z') ) {
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:143:7: ( '\"' (~ ( '\"' ) )* '\"' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:143:9: '\"' (~ ( '\"' ) )* '\"'
            {
            match('\"'); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:143:13: (~ ( '\"' ) )*
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

    // $ANTLR start "NESTED_ML_COMMENT"
    public final void mNESTED_ML_COMMENT() throws RecognitionException {
        try {
            int _type = NESTED_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:146:3: ( '[' (~ ( '[' | ']' ) )* ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:146:6: '[' (~ ( '[' | ']' ) )* ']'
            {
            match('['); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:146:10: (~ ( '[' | ']' ) )*
            loop3:
            do {
                int alt3=2;
                int LA3_0 = input.LA(1);

                if ( ((LA3_0 >= '\u0000' && LA3_0 <= 'Z')||LA3_0=='\\'||(LA3_0 >= '^' && LA3_0 <= '\uFFFF')) ) {
                    alt3=1;
                }


                switch (alt3) {
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
            	    break loop3;
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:149:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:149:9: ( ' ' | '\\t' | '\\r' | '\\n' )
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:156:6: ( . )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:156:8: .
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
        // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:8: ( T__20 | T__21 | T__22 | DEFAULT_INDICATOR | PHYLONET | TREE | UTREE | NETWORK | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | ID | QUOTE | NESTED_ML_COMMENT | WS | ELSE )
        int alt4=19;
        alt4 = dfa4.predict(input);
        switch (alt4) {
            case 1 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:10: T__20
                {
                mT__20(); 


                }
                break;
            case 2 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:16: T__21
                {
                mT__21(); 


                }
                break;
            case 3 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:22: T__22
                {
                mT__22(); 


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
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:124: NESTED_ML_COMMENT
                {
                mNESTED_ML_COMMENT(); 


                }
                break;
            case 18 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:142: WS
                {
                mWS(); 


                }
                break;
            case 19 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:1:145: ELSE
                {
                mELSE(); 


                }
                break;

        }

    }


    protected DFA4 dfa4 = new DFA4(this);
    static final String DFA4_eotS =
        "\5\uffff\4\26\1\20\2\26\1\uffff\2\20\6\uffff\1\26\1\uffff\3\26\1"+
        "\uffff\2\26\3\uffff\6\26\1\55\1\26\1\60\4\26\1\uffff\1\26\1\66\1"+
        "\uffff\1\26\1\70\1\26\1\72\1\26\1\uffff\1\26\1\uffff\1\26\1\uffff"+
        "\2\26\1\101\1\102\1\26\1\104\2\uffff\1\105\2\uffff";
    static final String DFA4_eofS =
        "\106\uffff";
    static final String DFA4_minS =
        "\1\0\4\uffff\1\110\1\122\1\124\1\105\1\116\1\105\1\116\1\uffff\2"+
        "\0\6\uffff\1\131\1\uffff\1\101\1\122\1\124\1\uffff\1\107\1\104\3"+
        "\uffff\1\114\1\105\1\116\1\105\1\127\1\111\1\60\1\117\1\60\1\123"+
        "\1\105\1\117\1\116\1\uffff\1\116\1\60\1\uffff\1\114\1\60\1\122\1"+
        "\60\1\105\1\uffff\1\101\1\uffff\1\113\1\uffff\2\124\2\60\1\105\1"+
        "\60\2\uffff\1\60\2\uffff";
    static final String DFA4_maxS =
        "\1\uffff\4\uffff\1\150\1\162\1\164\1\145\1\156\1\145\1\156\1\uffff"+
        "\2\uffff\6\uffff\1\171\1\uffff\1\145\1\162\1\164\1\uffff\1\147\1"+
        "\144\3\uffff\1\154\1\145\1\156\1\145\1\167\1\151\1\172\1\157\1\172"+
        "\1\163\1\145\1\157\1\156\1\uffff\1\156\1\172\1\uffff\1\154\1\172"+
        "\1\162\1\172\1\145\1\uffff\1\141\1\uffff\1\153\1\uffff\2\164\2\172"+
        "\1\145\1\172\2\uffff\1\172\2\uffff";
    static final String DFA4_acceptS =
        "\1\uffff\1\1\1\2\1\3\1\4\7\uffff\1\17\2\uffff\1\22\1\23\1\1\1\2"+
        "\1\3\1\4\1\uffff\1\17\3\uffff\1\11\2\uffff\1\20\1\21\1\22\15\uffff"+
        "\1\16\2\uffff\1\6\5\uffff\1\15\1\uffff\1\7\1\uffff\1\12\6\uffff"+
        "\1\10\1\5\1\uffff\1\14\1\13";
    static final String DFA4_specialS =
        "\1\2\14\uffff\1\1\1\0\67\uffff}>";
    static final String[] DFA4_transitionS = {
            "\11\20\2\17\2\20\1\17\22\20\1\17\1\20\1\15\1\11\6\20\1\4\1\20"+
            "\1\1\3\20\12\14\1\20\1\2\1\20\1\3\3\20\1\14\1\12\2\14\1\13\10"+
            "\14\1\10\1\14\1\5\3\14\1\6\1\7\5\14\1\16\5\20\1\14\1\12\2\14"+
            "\1\13\10\14\1\10\1\14\1\5\3\14\1\6\1\7\5\14\uff85\20",
            "",
            "",
            "",
            "",
            "\1\25\37\uffff\1\25",
            "\1\27\37\uffff\1\27",
            "\1\30\37\uffff\1\30",
            "\1\31\37\uffff\1\31",
            "\1\32\37\uffff\1\32",
            "\1\33\37\uffff\1\33",
            "\1\34\37\uffff\1\34",
            "",
            "\0\35",
            "\133\36\1\uffff\uffa4\36",
            "",
            "",
            "",
            "",
            "",
            "",
            "\1\40\37\uffff\1\40",
            "",
            "\1\42\3\uffff\1\41\33\uffff\1\42\3\uffff\1\41",
            "\1\43\37\uffff\1\43",
            "\1\44\37\uffff\1\44",
            "",
            "\1\45\37\uffff\1\45",
            "\1\46\37\uffff\1\46",
            "",
            "",
            "",
            "\1\47\37\uffff\1\47",
            "\1\50\37\uffff\1\50",
            "\1\51\37\uffff\1\51",
            "\1\52\37\uffff\1\52",
            "\1\53\37\uffff\1\53",
            "\1\54\37\uffff\1\54",
            "\12\26\7\uffff\32\26\6\uffff\32\26",
            "\1\56\37\uffff\1\56",
            "\12\26\7\uffff\22\26\1\57\7\26\6\uffff\22\26\1\57\7\26",
            "\1\61\37\uffff\1\61",
            "\1\62\37\uffff\1\62",
            "\1\63\37\uffff\1\63",
            "\1\64\37\uffff\1\64",
            "",
            "\1\65\37\uffff\1\65",
            "\12\26\7\uffff\32\26\6\uffff\32\26",
            "",
            "\1\67\37\uffff\1\67",
            "\12\26\7\uffff\32\26\6\uffff\32\26",
            "\1\71\37\uffff\1\71",
            "\12\26\7\uffff\32\26\6\uffff\32\26",
            "\1\73\37\uffff\1\73",
            "",
            "\1\74\37\uffff\1\74",
            "",
            "\1\75\37\uffff\1\75",
            "",
            "\1\76\37\uffff\1\76",
            "\1\77\37\uffff\1\77",
            "\12\26\7\uffff\22\26\1\100\7\26\6\uffff\22\26\1\100\7\26",
            "\12\26\7\uffff\32\26\6\uffff\32\26",
            "\1\103\37\uffff\1\103",
            "\12\26\7\uffff\32\26\6\uffff\32\26",
            "",
            "",
            "\12\26\7\uffff\32\26\6\uffff\32\26",
            "",
            ""
    };

    static final short[] DFA4_eot = DFA.unpackEncodedString(DFA4_eotS);
    static final short[] DFA4_eof = DFA.unpackEncodedString(DFA4_eofS);
    static final char[] DFA4_min = DFA.unpackEncodedStringToUnsignedChars(DFA4_minS);
    static final char[] DFA4_max = DFA.unpackEncodedStringToUnsignedChars(DFA4_maxS);
    static final short[] DFA4_accept = DFA.unpackEncodedString(DFA4_acceptS);
    static final short[] DFA4_special = DFA.unpackEncodedString(DFA4_specialS);
    static final short[][] DFA4_transition;

    static {
        int numStates = DFA4_transitionS.length;
        DFA4_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA4_transition[i] = DFA.unpackEncodedString(DFA4_transitionS[i]);
        }
    }

    class DFA4 extends DFA {

        public DFA4(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 4;
            this.eot = DFA4_eot;
            this.eof = DFA4_eof;
            this.min = DFA4_min;
            this.max = DFA4_max;
            this.accept = DFA4_accept;
            this.special = DFA4_special;
            this.transition = DFA4_transition;
        }
        public String getDescription() {
            return "1:1: Tokens : ( T__20 | T__21 | T__22 | DEFAULT_INDICATOR | PHYLONET | TREE | UTREE | NETWORK | START | BEGIN | TRANSLATE | NETWORKS | TREES | END | ID | QUOTE | NESTED_ML_COMMENT | WS | ELSE );";
        }
        public int specialStateTransition(int s, IntStream _input) throws NoViableAltException {
            IntStream input = _input;
        	int _s = s;
            switch ( s ) {
                    case 0 : 
                        int LA4_14 = input.LA(1);

                        s = -1;
                        if ( ((LA4_14 >= '\u0000' && LA4_14 <= 'Z')||(LA4_14 >= '\\' && LA4_14 <= '\uFFFF')) ) {s = 30;}

                        else s = 16;

                        if ( s>=0 ) return s;
                        break;

                    case 1 : 
                        int LA4_13 = input.LA(1);

                        s = -1;
                        if ( ((LA4_13 >= '\u0000' && LA4_13 <= '\uFFFF')) ) {s = 29;}

                        else s = 16;

                        if ( s>=0 ) return s;
                        break;

                    case 2 : 
                        int LA4_0 = input.LA(1);

                        s = -1;
                        if ( (LA4_0==',') ) {s = 1;}

                        else if ( (LA4_0==';') ) {s = 2;}

                        else if ( (LA4_0=='=') ) {s = 3;}

                        else if ( (LA4_0=='*') ) {s = 4;}

                        else if ( (LA4_0=='P'||LA4_0=='p') ) {s = 5;}

                        else if ( (LA4_0=='T'||LA4_0=='t') ) {s = 6;}

                        else if ( (LA4_0=='U'||LA4_0=='u') ) {s = 7;}

                        else if ( (LA4_0=='N'||LA4_0=='n') ) {s = 8;}

                        else if ( (LA4_0=='#') ) {s = 9;}

                        else if ( (LA4_0=='B'||LA4_0=='b') ) {s = 10;}

                        else if ( (LA4_0=='E'||LA4_0=='e') ) {s = 11;}

                        else if ( ((LA4_0 >= '0' && LA4_0 <= '9')||LA4_0=='A'||(LA4_0 >= 'C' && LA4_0 <= 'D')||(LA4_0 >= 'F' && LA4_0 <= 'M')||LA4_0=='O'||(LA4_0 >= 'Q' && LA4_0 <= 'S')||(LA4_0 >= 'V' && LA4_0 <= 'Z')||LA4_0=='a'||(LA4_0 >= 'c' && LA4_0 <= 'd')||(LA4_0 >= 'f' && LA4_0 <= 'm')||LA4_0=='o'||(LA4_0 >= 'q' && LA4_0 <= 's')||(LA4_0 >= 'v' && LA4_0 <= 'z')) ) {s = 12;}

                        else if ( (LA4_0=='\"') ) {s = 13;}

                        else if ( (LA4_0=='[') ) {s = 14;}

                        else if ( ((LA4_0 >= '\t' && LA4_0 <= '\n')||LA4_0=='\r'||LA4_0==' ') ) {s = 15;}

                        else if ( ((LA4_0 >= '\u0000' && LA4_0 <= '\b')||(LA4_0 >= '\u000B' && LA4_0 <= '\f')||(LA4_0 >= '\u000E' && LA4_0 <= '\u001F')||LA4_0=='!'||(LA4_0 >= '$' && LA4_0 <= ')')||LA4_0=='+'||(LA4_0 >= '-' && LA4_0 <= '/')||LA4_0==':'||LA4_0=='<'||(LA4_0 >= '>' && LA4_0 <= '@')||(LA4_0 >= '\\' && LA4_0 <= '`')||(LA4_0 >= '{' && LA4_0 <= '\uFFFF')) ) {s = 16;}

                        if ( s>=0 ) return s;
                        break;
            }
            NoViableAltException nvae =
                new NoViableAltException(getDescription(), 4, _s, input);
            error(nvae);
            throw nvae;
        }

    }
 

}