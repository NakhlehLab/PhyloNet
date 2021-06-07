// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g 2013-03-13 17:29:58

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.parsers.antlr.ast;


import org.antlr.runtime.*;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class RichNewick_1_1Lexer extends Lexer {
    public static final int EOF=-1;
    public static final int T__11=11;
    public static final int T__12=12;
    public static final int T__13=13;
    public static final int T__14=14;
    public static final int T__15=15;
    public static final int T__16=16;
    public static final int DECIMAL_NUMBER=4;
    public static final int NESTED_ML_COMMENT=5;
    public static final int QUOTED_TEXT=6;
    public static final int ROOTAGE_QUALIFIER=7;
    public static final int TREE_PROB=8;
    public static final int UNQUOTED_ALPHA_TEXT=9;
    public static final int WS=10;

    @Override
    public void emitErrorMessage(String msg)
    {
    // never want to dispaly errors to console, they may be recoverable
    }


    // delegates
    // delegators
    public Lexer[] getDelegates() {
        return new Lexer[] {};
    }

    public RichNewick_1_1Lexer() {} 
    public RichNewick_1_1Lexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public RichNewick_1_1Lexer(CharStream input, RecognizerSharedState state) {
        super(input,state);
    }
    public String getGrammarFileName() { return "D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g"; }

    // $ANTLR start "T__11"
    public final void mT__11() throws RecognitionException {
        try {
            int _type = T__11;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:13:7: ( '#' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:13:9: '#'
            {
            match('#'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__11"

    // $ANTLR start "T__12"
    public final void mT__12() throws RecognitionException {
        try {
            int _type = T__12;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:14:7: ( '(' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:14:9: '('
            {
            match('('); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__12"

    // $ANTLR start "T__13"
    public final void mT__13() throws RecognitionException {
        try {
            int _type = T__13;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:15:7: ( ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:15:9: ')'
            {
            match(')'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__13"

    // $ANTLR start "T__14"
    public final void mT__14() throws RecognitionException {
        try {
            int _type = T__14;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:16:7: ( ',' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:16:9: ','
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
    // $ANTLR end "T__14"

    // $ANTLR start "T__15"
    public final void mT__15() throws RecognitionException {
        try {
            int _type = T__15;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:17:7: ( ':' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:17:9: ':'
            {
            match(':'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "T__15"

    // $ANTLR start "T__16"
    public final void mT__16() throws RecognitionException {
        try {
            int _type = T__16;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:18:7: ( ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:18:9: ';'
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
    // $ANTLR end "T__16"

    // $ANTLR start "ROOTAGE_QUALIFIER"
    public final void mROOTAGE_QUALIFIER() throws RecognitionException {
        try {
            int _type = ROOTAGE_QUALIFIER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:117:2: ( '[' '&' ( 'r' | 'R' | 'u' | 'U' ) ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:117:4: '[' '&' ( 'r' | 'R' | 'u' | 'U' ) ']'
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

    // $ANTLR start "TREE_PROB"
    public final void mTREE_PROB() throws RecognitionException {
        try {
            int _type = TREE_PROB;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:2: ( '[' '&' ( 'W' | 'w' ) ( ' ' )+ ( ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ ) | ( '0' .. '9' )+ ) ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:4: '[' '&' ( 'W' | 'w' ) ( ' ' )+ ( ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ ) | ( '0' .. '9' )+ ) ']'
            {
            match('['); 

            match('&'); 

            if ( input.LA(1)=='W'||input.LA(1)=='w' ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:22: ( ' ' )+
            int cnt1=0;
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==' ') ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:22: ' '
            	    {
            	    match(' '); 

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


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:27: ( ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ ) | ( '0' .. '9' )+ )
            int alt5=2;
            alt5 = dfa5.predict(input);
            switch (alt5) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:28: ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ )
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:28: ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ )
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:29: ( '0' .. '9' )* '.' ( '0' .. '9' )+
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:29: ( '0' .. '9' )*
                    loop2:
                    do {
                        int alt2=2;
                        int LA2_0 = input.LA(1);

                        if ( ((LA2_0 >= '0' && LA2_0 <= '9')) ) {
                            alt2=1;
                        }


                        switch (alt2) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
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


                    match('.'); 

                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:45: ( '0' .. '9' )+
                    int cnt3=0;
                    loop3:
                    do {
                        int alt3=2;
                        int LA3_0 = input.LA(1);

                        if ( ((LA3_0 >= '0' && LA3_0 <= '9')) ) {
                            alt3=1;
                        }


                        switch (alt3) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
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
                    	    if ( cnt3 >= 1 ) break loop3;
                                EarlyExitException eee =
                                    new EarlyExitException(3, input);
                                throw eee;
                        }
                        cnt3++;
                    } while (true);


                    }


                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:58: ( '0' .. '9' )+
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:121:58: ( '0' .. '9' )+
                    int cnt4=0;
                    loop4:
                    do {
                        int alt4=2;
                        int LA4_0 = input.LA(1);

                        if ( ((LA4_0 >= '0' && LA4_0 <= '9')) ) {
                            alt4=1;
                        }


                        switch (alt4) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
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
                    break;

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
    // $ANTLR end "TREE_PROB"

    // $ANTLR start "DECIMAL_NUMBER"
    public final void mDECIMAL_NUMBER() throws RecognitionException {
        try {
            int _type = DECIMAL_NUMBER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:16: ( ( ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ ) | ( '0' .. '9' )+ ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:18: ( ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ ) | ( '0' .. '9' )+ )
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:18: ( ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ ) | ( '0' .. '9' )+ )
            int alt9=2;
            alt9 = dfa9.predict(input);
            switch (alt9) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:19: ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ )
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:19: ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ )
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:20: ( '0' .. '9' )* '.' ( '0' .. '9' )+
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:20: ( '0' .. '9' )*
                    loop6:
                    do {
                        int alt6=2;
                        int LA6_0 = input.LA(1);

                        if ( ((LA6_0 >= '0' && LA6_0 <= '9')) ) {
                            alt6=1;
                        }


                        switch (alt6) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
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
                    	    break loop6;
                        }
                    } while (true);


                    match('.'); 

                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:36: ( '0' .. '9' )+
                    int cnt7=0;
                    loop7:
                    do {
                        int alt7=2;
                        int LA7_0 = input.LA(1);

                        if ( ((LA7_0 >= '0' && LA7_0 <= '9')) ) {
                            alt7=1;
                        }


                        switch (alt7) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
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
                    	    if ( cnt7 >= 1 ) break loop7;
                                EarlyExitException eee =
                                    new EarlyExitException(7, input);
                                throw eee;
                        }
                        cnt7++;
                    } while (true);


                    }


                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:49: ( '0' .. '9' )+
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:124:49: ( '0' .. '9' )+
                    int cnt8=0;
                    loop8:
                    do {
                        int alt8=2;
                        int LA8_0 = input.LA(1);

                        if ( ((LA8_0 >= '0' && LA8_0 <= '9')) ) {
                            alt8=1;
                        }


                        switch (alt8) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:
                    	    {
                    	    if ( (input.LA(1) >= '0' && input.LA(1) <= '9') ) {
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
                    break;

            }


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "DECIMAL_NUMBER"

    // $ANTLR start "UNQUOTED_ALPHA_TEXT"
    public final void mUNQUOTED_ALPHA_TEXT() throws RecognitionException {
        try {
            int _type = UNQUOTED_ALPHA_TEXT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:127:2: (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | ' ' | '.' | ( '0' .. '9' ) | '\\t' | '\\r' | '\\n' ) (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | '.' | ( '0' .. '9' ) ) )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:127:4: ~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | ' ' | '.' | ( '0' .. '9' ) | '\\t' | '\\r' | '\\n' ) (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | '.' | ( '0' .. '9' ) ) )*
            {
            if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '\b')||(input.LA(1) >= '\u000B' && input.LA(1) <= '\f')||(input.LA(1) >= '\u000E' && input.LA(1) <= '\u001F')||(input.LA(1) >= '!' && input.LA(1) <= '\"')||(input.LA(1) >= '$' && input.LA(1) <= '&')||(input.LA(1) >= '*' && input.LA(1) <= '+')||input.LA(1)=='-'||input.LA(1)=='/'||(input.LA(1) >= '<' && input.LA(1) <= 'Z')||input.LA(1)=='\\'||(input.LA(1) >= '^' && input.LA(1) <= '\uFFFF') ) {
                input.consume();
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;
            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:128:32: (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | '.' | ( '0' .. '9' ) ) )*
            loop10:
            do {
                int alt10=2;
                int LA10_0 = input.LA(1);

                if ( ((LA10_0 >= '\u0000' && LA10_0 <= '\"')||(LA10_0 >= '$' && LA10_0 <= '&')||(LA10_0 >= '*' && LA10_0 <= '+')||LA10_0=='-'||LA10_0=='/'||(LA10_0 >= '<' && LA10_0 <= 'Z')||LA10_0=='\\'||(LA10_0 >= '^' && LA10_0 <= '\uFFFF')) ) {
                    alt10=1;
                }


                switch (alt10) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '\"')||(input.LA(1) >= '$' && input.LA(1) <= '&')||(input.LA(1) >= '*' && input.LA(1) <= '+')||input.LA(1)=='-'||input.LA(1)=='/'||(input.LA(1) >= '<' && input.LA(1) <= 'Z')||input.LA(1)=='\\'||(input.LA(1) >= '^' && input.LA(1) <= '\uFFFF') ) {
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
            	    break loop10;
                }
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "UNQUOTED_ALPHA_TEXT"

    // $ANTLR start "QUOTED_TEXT"
    public final void mQUOTED_TEXT() throws RecognitionException {
        try {
            int _type = QUOTED_TEXT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:131:2: ( '\\'' (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )* '\\'' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:131:4: '\\'' (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )* '\\''
            {
            match('\''); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:131:9: (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )*
            loop11:
            do {
                int alt11=3;
                int LA11_0 = input.LA(1);

                if ( (LA11_0=='\'') ) {
                    int LA11_1 = input.LA(2);

                    if ( (LA11_1=='\'') ) {
                        alt11=2;
                    }


                }
                else if ( ((LA11_0 >= '\u0000' && LA11_0 <= '\t')||(LA11_0 >= '\u000B' && LA11_0 <= '\f')||(LA11_0 >= '\u000E' && LA11_0 <= '&')||(LA11_0 >= '(' && LA11_0 <= '\uFFFF')) ) {
                    alt11=1;
                }


                switch (alt11) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:131:10: ~ ( '\\n' | '\\r' | '\\'' )
            	    {
            	    if ( (input.LA(1) >= '\u0000' && input.LA(1) <= '\t')||(input.LA(1) >= '\u000B' && input.LA(1) <= '\f')||(input.LA(1) >= '\u000E' && input.LA(1) <= '&')||(input.LA(1) >= '(' && input.LA(1) <= '\uFFFF') ) {
            	        input.consume();
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;
            	    }


            	    }
            	    break;
            	case 2 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:131:34: ( '\\'\\'' )
            	    {
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:131:34: ( '\\'\\'' )
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:131:35: '\\'\\''
            	    {
            	    match("''"); 



            	    }


            	    }
            	    break;

            	default :
            	    break loop11;
                }
            } while (true);


            match('\''); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        	// do for sure before leaving
        }
    }
    // $ANTLR end "QUOTED_TEXT"

    // $ANTLR start "NESTED_ML_COMMENT"
    public final void mNESTED_ML_COMMENT() throws RecognitionException {
        try {
            int _type = NESTED_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:134:3: ( '[' (~ ( '[' | ']' ) | NESTED_ML_COMMENT )* ']' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:134:6: '[' (~ ( '[' | ']' ) | NESTED_ML_COMMENT )* ']'
            {
            match('['); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:134:10: (~ ( '[' | ']' ) | NESTED_ML_COMMENT )*
            loop12:
            do {
                int alt12=3;
                int LA12_0 = input.LA(1);

                if ( ((LA12_0 >= '\u0000' && LA12_0 <= 'Z')||LA12_0=='\\'||(LA12_0 >= '^' && LA12_0 <= '\uFFFF')) ) {
                    alt12=1;
                }
                else if ( (LA12_0=='[') ) {
                    alt12=2;
                }


                switch (alt12) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:134:11: ~ ( '[' | ']' )
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
            	case 2 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:134:26: NESTED_ML_COMMENT
            	    {
            	    mNESTED_ML_COMMENT(); 


            	    }
            	    break;

            	default :
            	    break loop12;
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
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:136:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:136:9: ( ' ' | '\\t' | '\\r' | '\\n' )
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

    public void mTokens() throws RecognitionException {
        // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:8: ( T__11 | T__12 | T__13 | T__14 | T__15 | T__16 | ROOTAGE_QUALIFIER | TREE_PROB | DECIMAL_NUMBER | UNQUOTED_ALPHA_TEXT | QUOTED_TEXT | NESTED_ML_COMMENT | WS )
        int alt13=13;
        alt13 = dfa13.predict(input);
        switch (alt13) {
            case 1 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:10: T__11
                {
                mT__11(); 


                }
                break;
            case 2 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:16: T__12
                {
                mT__12(); 


                }
                break;
            case 3 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:22: T__13
                {
                mT__13(); 


                }
                break;
            case 4 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:28: T__14
                {
                mT__14(); 


                }
                break;
            case 5 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:34: T__15
                {
                mT__15(); 


                }
                break;
            case 6 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:40: T__16
                {
                mT__16(); 


                }
                break;
            case 7 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:46: ROOTAGE_QUALIFIER
                {
                mROOTAGE_QUALIFIER(); 


                }
                break;
            case 8 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:64: TREE_PROB
                {
                mTREE_PROB(); 


                }
                break;
            case 9 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:74: DECIMAL_NUMBER
                {
                mDECIMAL_NUMBER(); 


                }
                break;
            case 10 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:89: UNQUOTED_ALPHA_TEXT
                {
                mUNQUOTED_ALPHA_TEXT(); 


                }
                break;
            case 11 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:109: QUOTED_TEXT
                {
                mQUOTED_TEXT(); 


                }
                break;
            case 12 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:121: NESTED_ML_COMMENT
                {
                mNESTED_ML_COMMENT(); 


                }
                break;
            case 13 :
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:1:139: WS
                {
                mWS(); 


                }
                break;

        }

    }


    protected DFA5 dfa5 = new DFA5(this);
    protected DFA9 dfa9 = new DFA9(this);
    protected DFA13 dfa13 = new DFA13(this);
    static final String DFA5_eotS =
        "\4\uffff";
    static final String DFA5_eofS =
        "\4\uffff";
    static final String DFA5_minS =
        "\2\56\2\uffff";
    static final String DFA5_maxS =
        "\1\71\1\135\2\uffff";
    static final String DFA5_acceptS =
        "\2\uffff\1\1\1\2";
    static final String DFA5_specialS =
        "\4\uffff}>";
    static final String[] DFA5_transitionS = {
            "\1\2\1\uffff\12\1",
            "\1\2\1\uffff\12\1\43\uffff\1\3",
            "",
            ""
    };

    static final short[] DFA5_eot = DFA.unpackEncodedString(DFA5_eotS);
    static final short[] DFA5_eof = DFA.unpackEncodedString(DFA5_eofS);
    static final char[] DFA5_min = DFA.unpackEncodedStringToUnsignedChars(DFA5_minS);
    static final char[] DFA5_max = DFA.unpackEncodedStringToUnsignedChars(DFA5_maxS);
    static final short[] DFA5_accept = DFA.unpackEncodedString(DFA5_acceptS);
    static final short[] DFA5_special = DFA.unpackEncodedString(DFA5_specialS);
    static final short[][] DFA5_transition;

    static {
        int numStates = DFA5_transitionS.length;
        DFA5_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA5_transition[i] = DFA.unpackEncodedString(DFA5_transitionS[i]);
        }
    }

    class DFA5 extends DFA {

        public DFA5(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 5;
            this.eot = DFA5_eot;
            this.eof = DFA5_eof;
            this.min = DFA5_min;
            this.max = DFA5_max;
            this.accept = DFA5_accept;
            this.special = DFA5_special;
            this.transition = DFA5_transition;
        }
        public String getDescription() {
            return "121:27: ( ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ ) | ( '0' .. '9' )+ )";
        }
    }
    static final String DFA9_eotS =
        "\1\uffff\1\3\2\uffff";
    static final String DFA9_eofS =
        "\4\uffff";
    static final String DFA9_minS =
        "\2\56\2\uffff";
    static final String DFA9_maxS =
        "\2\71\2\uffff";
    static final String DFA9_acceptS =
        "\2\uffff\1\1\1\2";
    static final String DFA9_specialS =
        "\4\uffff}>";
    static final String[] DFA9_transitionS = {
            "\1\2\1\uffff\12\1",
            "\1\2\1\uffff\12\1",
            "",
            ""
    };

    static final short[] DFA9_eot = DFA.unpackEncodedString(DFA9_eotS);
    static final short[] DFA9_eof = DFA.unpackEncodedString(DFA9_eofS);
    static final char[] DFA9_min = DFA.unpackEncodedStringToUnsignedChars(DFA9_minS);
    static final char[] DFA9_max = DFA.unpackEncodedStringToUnsignedChars(DFA9_maxS);
    static final short[] DFA9_accept = DFA.unpackEncodedString(DFA9_acceptS);
    static final short[] DFA9_special = DFA.unpackEncodedString(DFA9_specialS);
    static final short[][] DFA9_transition;

    static {
        int numStates = DFA9_transitionS.length;
        DFA9_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA9_transition[i] = DFA.unpackEncodedString(DFA9_transitionS[i]);
        }
    }

    class DFA9 extends DFA {

        public DFA9(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 9;
            this.eot = DFA9_eot;
            this.eof = DFA9_eof;
            this.min = DFA9_min;
            this.max = DFA9_max;
            this.accept = DFA9_accept;
            this.special = DFA9_special;
            this.transition = DFA9_transition;
        }
        public String getDescription() {
            return "124:18: ( ( ( '0' .. '9' )* '.' ( '0' .. '9' )+ ) | ( '0' .. '9' )+ )";
        }
    }
    static final String DFA13_eotS =
        "\30\uffff";
    static final String DFA13_eofS =
        "\30\uffff";
    static final String DFA13_minS =
        "\1\0\6\uffff\1\0\4\uffff\1\0\1\uffff\2\0\1\uffff\1\0\1\uffff\2\0"+
        "\1\uffff\1\0\1\uffff";
    static final String DFA13_maxS =
        "\1\uffff\6\uffff\1\uffff\4\uffff\1\uffff\1\uffff\2\uffff\1\uffff"+
        "\1\uffff\1\uffff\2\uffff\1\uffff\1\uffff\1\uffff";
    static final String DFA13_acceptS =
        "\1\uffff\1\1\1\2\1\3\1\4\1\5\1\6\1\uffff\1\11\1\12\1\13\1\15\1\uffff"+
        "\1\14\2\uffff\1\7\1\uffff\1\7\2\uffff\1\10\1\uffff\1\10";
    static final String DFA13_specialS =
        "\1\2\6\uffff\1\6\4\uffff\1\5\1\uffff\1\10\1\0\1\uffff\1\3\1\uffff"+
        "\1\1\1\4\1\uffff\1\7\1\uffff}>";
    static final String[] DFA13_transitionS = {
            "\11\11\2\13\2\11\1\13\22\11\1\13\2\11\1\1\3\11\1\12\1\2\1\3"+
            "\2\11\1\4\1\11\1\10\1\11\12\10\1\5\1\6\37\11\1\7\1\11\1\uffff"+
            "\uffa2\11",
            "",
            "",
            "",
            "",
            "",
            "",
            "\46\15\1\14\uffd9\15",
            "",
            "",
            "",
            "",
            "\122\15\1\16\2\15\1\16\1\15\1\17\32\15\1\16\2\15\1\16\1\15"+
            "\1\17\uff88\15",
            "",
            "\135\15\1\20\uffa2\15",
            "\40\15\1\21\uffdf\15",
            "",
            "\40\15\1\21\15\15\1\24\1\15\12\23\uffc6\15",
            "",
            "\56\15\1\24\1\15\12\23\43\15\1\25\uffa2\15",
            "\60\15\12\26\uffc6\15",
            "",
            "\60\15\12\26\43\15\1\25\uffa2\15",
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
            return "1:1: Tokens : ( T__11 | T__12 | T__13 | T__14 | T__15 | T__16 | ROOTAGE_QUALIFIER | TREE_PROB | DECIMAL_NUMBER | UNQUOTED_ALPHA_TEXT | QUOTED_TEXT | NESTED_ML_COMMENT | WS );";
        }
        public int specialStateTransition(int s, IntStream _input) throws NoViableAltException {
            IntStream input = _input;
        	int _s = s;
            switch ( s ) {
                    case 0 : 
                        int LA13_15 = input.LA(1);

                        s = -1;
                        if ( (LA13_15==' ') ) {s = 17;}

                        else if ( ((LA13_15 >= '\u0000' && LA13_15 <= '\u001F')||(LA13_15 >= '!' && LA13_15 <= '\uFFFF')) ) {s = 13;}

                        if ( s>=0 ) return s;
                        break;

                    case 1 : 
                        int LA13_19 = input.LA(1);

                        s = -1;
                        if ( (LA13_19=='.') ) {s = 20;}

                        else if ( ((LA13_19 >= '0' && LA13_19 <= '9')) ) {s = 19;}

                        else if ( (LA13_19==']') ) {s = 21;}

                        else if ( ((LA13_19 >= '\u0000' && LA13_19 <= '-')||LA13_19=='/'||(LA13_19 >= ':' && LA13_19 <= '\\')||(LA13_19 >= '^' && LA13_19 <= '\uFFFF')) ) {s = 13;}

                        if ( s>=0 ) return s;
                        break;

                    case 2 : 
                        int LA13_0 = input.LA(1);

                        s = -1;
                        if ( (LA13_0=='#') ) {s = 1;}

                        else if ( (LA13_0=='(') ) {s = 2;}

                        else if ( (LA13_0==')') ) {s = 3;}

                        else if ( (LA13_0==',') ) {s = 4;}

                        else if ( (LA13_0==':') ) {s = 5;}

                        else if ( (LA13_0==';') ) {s = 6;}

                        else if ( (LA13_0=='[') ) {s = 7;}

                        else if ( (LA13_0=='.'||(LA13_0 >= '0' && LA13_0 <= '9')) ) {s = 8;}

                        else if ( ((LA13_0 >= '\u0000' && LA13_0 <= '\b')||(LA13_0 >= '\u000B' && LA13_0 <= '\f')||(LA13_0 >= '\u000E' && LA13_0 <= '\u001F')||(LA13_0 >= '!' && LA13_0 <= '\"')||(LA13_0 >= '$' && LA13_0 <= '&')||(LA13_0 >= '*' && LA13_0 <= '+')||LA13_0=='-'||LA13_0=='/'||(LA13_0 >= '<' && LA13_0 <= 'Z')||LA13_0=='\\'||(LA13_0 >= '^' && LA13_0 <= '\uFFFF')) ) {s = 9;}

                        else if ( (LA13_0=='\'') ) {s = 10;}

                        else if ( ((LA13_0 >= '\t' && LA13_0 <= '\n')||LA13_0=='\r'||LA13_0==' ') ) {s = 11;}

                        if ( s>=0 ) return s;
                        break;

                    case 3 : 
                        int LA13_17 = input.LA(1);

                        s = -1;
                        if ( ((LA13_17 >= '0' && LA13_17 <= '9')) ) {s = 19;}

                        else if ( (LA13_17=='.') ) {s = 20;}

                        else if ( (LA13_17==' ') ) {s = 17;}

                        else if ( ((LA13_17 >= '\u0000' && LA13_17 <= '\u001F')||(LA13_17 >= '!' && LA13_17 <= '-')||LA13_17=='/'||(LA13_17 >= ':' && LA13_17 <= '\uFFFF')) ) {s = 13;}

                        if ( s>=0 ) return s;
                        break;

                    case 4 : 
                        int LA13_20 = input.LA(1);

                        s = -1;
                        if ( ((LA13_20 >= '0' && LA13_20 <= '9')) ) {s = 22;}

                        else if ( ((LA13_20 >= '\u0000' && LA13_20 <= '/')||(LA13_20 >= ':' && LA13_20 <= '\uFFFF')) ) {s = 13;}

                        if ( s>=0 ) return s;
                        break;

                    case 5 : 
                        int LA13_12 = input.LA(1);

                        s = -1;
                        if ( (LA13_12=='R'||LA13_12=='U'||LA13_12=='r'||LA13_12=='u') ) {s = 14;}

                        else if ( (LA13_12=='W'||LA13_12=='w') ) {s = 15;}

                        else if ( ((LA13_12 >= '\u0000' && LA13_12 <= 'Q')||(LA13_12 >= 'S' && LA13_12 <= 'T')||LA13_12=='V'||(LA13_12 >= 'X' && LA13_12 <= 'q')||(LA13_12 >= 's' && LA13_12 <= 't')||LA13_12=='v'||(LA13_12 >= 'x' && LA13_12 <= '\uFFFF')) ) {s = 13;}

                        if ( s>=0 ) return s;
                        break;

                    case 6 : 
                        int LA13_7 = input.LA(1);

                        s = -1;
                        if ( (LA13_7=='&') ) {s = 12;}

                        else if ( ((LA13_7 >= '\u0000' && LA13_7 <= '%')||(LA13_7 >= '\'' && LA13_7 <= '\uFFFF')) ) {s = 13;}

                        if ( s>=0 ) return s;
                        break;

                    case 7 : 
                        int LA13_22 = input.LA(1);

                        s = -1;
                        if ( (LA13_22==']') ) {s = 21;}

                        else if ( ((LA13_22 >= '0' && LA13_22 <= '9')) ) {s = 22;}

                        else if ( ((LA13_22 >= '\u0000' && LA13_22 <= '/')||(LA13_22 >= ':' && LA13_22 <= '\\')||(LA13_22 >= '^' && LA13_22 <= '\uFFFF')) ) {s = 13;}

                        if ( s>=0 ) return s;
                        break;

                    case 8 : 
                        int LA13_14 = input.LA(1);

                        s = -1;
                        if ( (LA13_14==']') ) {s = 16;}

                        else if ( ((LA13_14 >= '\u0000' && LA13_14 <= '\\')||(LA13_14 >= '^' && LA13_14 <= '\uFFFF')) ) {s = 13;}

                        if ( s>=0 ) return s;
                        break;
            }
            NoViableAltException nvae =
                new NoViableAltException(getDescription(), 13, _s, input);
            error(nvae);
            throw nvae;
        }

    }
 

}