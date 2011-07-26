// $ANTLR 3.3 Nov 30, 2010 12:45:30 C:\\Users\\Matt\\Desktop\\ExtendedNewick.g 2011-07-21 15:22:27

package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

public class ExtendedNewickLexer extends Lexer {
    public static final int EOF=-1;
    public static final int T__9=9;
    public static final int T__10=10;
    public static final int T__11=11;
    public static final int T__12=12;
    public static final int T__13=13;
    public static final int T__14=14;
    public static final int UNQUOTED_ALPHA_TEXT=4;
    public static final int DECIMAL_NUMBER=5;
    public static final int QUOTED_TEXT=6;
    public static final int NESTED_ML_COMMENT=7;
    public static final int WS=8;

    // delegates
    // delegators

    public ExtendedNewickLexer() {;} 
    public ExtendedNewickLexer(CharStream input) {
        this(input, new RecognizerSharedState());
    }
    public ExtendedNewickLexer(CharStream input, RecognizerSharedState state) {
        super(input,state);

    }
    public String getGrammarFileName() { return "C:\\Users\\Matt\\Desktop\\ExtendedNewick.g"; }

    // $ANTLR start "T__9"
    public final void mT__9() throws RecognitionException {
        try {
            int _type = T__9;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:7:6: ( ';' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:7:8: ';'
            {
            match(';'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__9"

    // $ANTLR start "T__10"
    public final void mT__10() throws RecognitionException {
        try {
            int _type = T__10;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:8:7: ( '(' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:8:9: '('
            {
            match('('); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__10"

    // $ANTLR start "T__11"
    public final void mT__11() throws RecognitionException {
        try {
            int _type = T__11;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:9:7: ( ',' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:9:9: ','
            {
            match(','); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__11"

    // $ANTLR start "T__12"
    public final void mT__12() throws RecognitionException {
        try {
            int _type = T__12;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:10:7: ( ')' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:10:9: ')'
            {
            match(')'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__12"

    // $ANTLR start "T__13"
    public final void mT__13() throws RecognitionException {
        try {
            int _type = T__13;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:11:7: ( ':' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:11:9: ':'
            {
            match(':'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__13"

    // $ANTLR start "T__14"
    public final void mT__14() throws RecognitionException {
        try {
            int _type = T__14;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:12:7: ( '#' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:12:9: '#'
            {
            match('#'); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "T__14"

    // $ANTLR start "NESTED_ML_COMMENT"
    public final void mNESTED_ML_COMMENT() throws RecognitionException {
        try {
            int _type = NESTED_ML_COMMENT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:81:3: ( '[' (~ ( '[' | ']' ) | NESTED_ML_COMMENT )* ']' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:81:6: '[' (~ ( '[' | ']' ) | NESTED_ML_COMMENT )* ']'
            {
            match('['); 
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:81:10: (~ ( '[' | ']' ) | NESTED_ML_COMMENT )*
            loop1:
            do {
                int alt1=3;
                int LA1_0 = input.LA(1);

                if ( ((LA1_0>='\u0000' && LA1_0<='Z')||LA1_0=='\\'||(LA1_0>='^' && LA1_0<='\uFFFF')) ) {
                    alt1=1;
                }
                else if ( (LA1_0=='[') ) {
                    alt1=2;
                }


                switch (alt1) {
            	case 1 :
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:81:11: ~ ( '[' | ']' )
            	    {
            	    if ( (input.LA(1)>='\u0000' && input.LA(1)<='Z')||input.LA(1)=='\\'||(input.LA(1)>='^' && input.LA(1)<='\uFFFF') ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;
            	case 2 :
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:81:26: NESTED_ML_COMMENT
            	    {
            	    mNESTED_ML_COMMENT(); 

            	    }
            	    break;

            	default :
            	    break loop1;
                }
            } while (true);

            match(']'); 
            _channel=HIDDEN;

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "NESTED_ML_COMMENT"

    // $ANTLR start "WS"
    public final void mWS() throws RecognitionException {
        try {
            int _type = WS;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:83:5: ( ( ' ' | '\\t' | '\\r' | '\\n' ) )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:83:9: ( ' ' | '\\t' | '\\r' | '\\n' )
            {
            if ( (input.LA(1)>='\t' && input.LA(1)<='\n')||input.LA(1)=='\r'||input.LA(1)==' ' ) {
                input.consume();

            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                recover(mse);
                throw mse;}

            _channel=HIDDEN;

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "WS"

    // $ANTLR start "DECIMAL_NUMBER"
    public final void mDECIMAL_NUMBER() throws RecognitionException {
        try {
            int _type = DECIMAL_NUMBER;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:89:16: ( ( ( '0' .. '9' ) | '.' )+ )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:89:18: ( ( '0' .. '9' ) | '.' )+
            {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:89:18: ( ( '0' .. '9' ) | '.' )+
            int cnt2=0;
            loop2:
            do {
                int alt2=3;
                int LA2_0 = input.LA(1);

                if ( ((LA2_0>='0' && LA2_0<='9')) ) {
                    alt2=1;
                }
                else if ( (LA2_0=='.') ) {
                    alt2=2;
                }


                switch (alt2) {
            	case 1 :
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:89:19: ( '0' .. '9' )
            	    {
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:89:19: ( '0' .. '9' )
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:89:20: '0' .. '9'
            	    {
            	    matchRange('0','9'); 

            	    }


            	    }
            	    break;
            	case 2 :
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:89:32: '.'
            	    {
            	    match('.'); 

            	    }
            	    break;

            	default :
            	    if ( cnt2 >= 1 ) break loop2;
                        EarlyExitException eee =
                            new EarlyExitException(2, input);
                        throw eee;
                }
                cnt2++;
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "DECIMAL_NUMBER"

    // $ANTLR start "UNQUOTED_ALPHA_TEXT"
    public final void mUNQUOTED_ALPHA_TEXT() throws RecognitionException {
        try {
            int _type = UNQUOTED_ALPHA_TEXT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:92:2: ( (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | ' ' | '.' | ( '0' .. '9' ) | '\\t' | '\\r' | '\\n' ) )* )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:92:4: (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | ' ' | '.' | ( '0' .. '9' ) | '\\t' | '\\r' | '\\n' ) )*
            {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:92:4: (~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | ' ' | '.' | ( '0' .. '9' ) | '\\t' | '\\r' | '\\n' ) )*
            loop3:
            do {
                int alt3=2;
                int LA3_0 = input.LA(1);

                if ( ((LA3_0>='\u0000' && LA3_0<='\b')||(LA3_0>='\u000B' && LA3_0<='\f')||(LA3_0>='\u000E' && LA3_0<='\u001F')||(LA3_0>='!' && LA3_0<='\"')||(LA3_0>='$' && LA3_0<='&')||(LA3_0>='*' && LA3_0<='+')||LA3_0=='-'||LA3_0=='/'||(LA3_0>='<' && LA3_0<='Z')||LA3_0=='\\'||(LA3_0>='^' && LA3_0<='\uFFFF')) ) {
                    alt3=1;
                }


                switch (alt3) {
            	case 1 :
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:92:4: ~ ( '(' | ')' | '[' | ']' | ':' | ';' | '#' | '\\'' | ',' | ' ' | '.' | ( '0' .. '9' ) | '\\t' | '\\r' | '\\n' )
            	    {
            	    if ( (input.LA(1)>='\u0000' && input.LA(1)<='\b')||(input.LA(1)>='\u000B' && input.LA(1)<='\f')||(input.LA(1)>='\u000E' && input.LA(1)<='\u001F')||(input.LA(1)>='!' && input.LA(1)<='\"')||(input.LA(1)>='$' && input.LA(1)<='&')||(input.LA(1)>='*' && input.LA(1)<='+')||input.LA(1)=='-'||input.LA(1)=='/'||(input.LA(1)>='<' && input.LA(1)<='Z')||input.LA(1)=='\\'||(input.LA(1)>='^' && input.LA(1)<='\uFFFF') ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;

            	default :
            	    break loop3;
                }
            } while (true);


            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "UNQUOTED_ALPHA_TEXT"

    // $ANTLR start "QUOTED_TEXT"
    public final void mQUOTED_TEXT() throws RecognitionException {
        try {
            int _type = QUOTED_TEXT;
            int _channel = DEFAULT_TOKEN_CHANNEL;
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:96:2: ( '\\'' (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )* '\\'' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:96:4: '\\'' (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )* '\\''
            {
            match('\''); 
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:96:9: (~ ( '\\n' | '\\r' | '\\'' ) | ( '\\'\\'' ) )*
            loop4:
            do {
                int alt4=3;
                int LA4_0 = input.LA(1);

                if ( (LA4_0=='\'') ) {
                    int LA4_1 = input.LA(2);

                    if ( (LA4_1=='\'') ) {
                        alt4=2;
                    }


                }
                else if ( ((LA4_0>='\u0000' && LA4_0<='\t')||(LA4_0>='\u000B' && LA4_0<='\f')||(LA4_0>='\u000E' && LA4_0<='&')||(LA4_0>='(' && LA4_0<='\uFFFF')) ) {
                    alt4=1;
                }


                switch (alt4) {
            	case 1 :
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:96:10: ~ ( '\\n' | '\\r' | '\\'' )
            	    {
            	    if ( (input.LA(1)>='\u0000' && input.LA(1)<='\t')||(input.LA(1)>='\u000B' && input.LA(1)<='\f')||(input.LA(1)>='\u000E' && input.LA(1)<='&')||(input.LA(1)>='(' && input.LA(1)<='\uFFFF') ) {
            	        input.consume();

            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        recover(mse);
            	        throw mse;}


            	    }
            	    break;
            	case 2 :
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:96:34: ( '\\'\\'' )
            	    {
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:96:34: ( '\\'\\'' )
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:96:35: '\\'\\''
            	    {
            	    match("''"); 


            	    }


            	    }
            	    break;

            	default :
            	    break loop4;
                }
            } while (true);

            match('\''); 

            }

            state.type = _type;
            state.channel = _channel;
        }
        finally {
        }
    }
    // $ANTLR end "QUOTED_TEXT"

    public void mTokens() throws RecognitionException {
        // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:8: ( T__9 | T__10 | T__11 | T__12 | T__13 | T__14 | NESTED_ML_COMMENT | WS | DECIMAL_NUMBER | UNQUOTED_ALPHA_TEXT | QUOTED_TEXT )
        int alt5=11;
        alt5 = dfa5.predict(input);
        switch (alt5) {
            case 1 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:10: T__9
                {
                mT__9(); 

                }
                break;
            case 2 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:15: T__10
                {
                mT__10(); 

                }
                break;
            case 3 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:21: T__11
                {
                mT__11(); 

                }
                break;
            case 4 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:27: T__12
                {
                mT__12(); 

                }
                break;
            case 5 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:33: T__13
                {
                mT__13(); 

                }
                break;
            case 6 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:39: T__14
                {
                mT__14(); 

                }
                break;
            case 7 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:45: NESTED_ML_COMMENT
                {
                mNESTED_ML_COMMENT(); 

                }
                break;
            case 8 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:63: WS
                {
                mWS(); 

                }
                break;
            case 9 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:66: DECIMAL_NUMBER
                {
                mDECIMAL_NUMBER(); 

                }
                break;
            case 10 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:81: UNQUOTED_ALPHA_TEXT
                {
                mUNQUOTED_ALPHA_TEXT(); 

                }
                break;
            case 11 :
                // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:1:101: QUOTED_TEXT
                {
                mQUOTED_TEXT(); 

                }
                break;

        }

    }


    protected DFA5 dfa5 = new DFA5(this);
    static final String DFA5_eotS =
        "\1\12\13\uffff";
    static final String DFA5_eofS =
        "\14\uffff";
    static final String DFA5_minS =
        "\1\11\13\uffff";
    static final String DFA5_maxS =
        "\1\133\13\uffff";
    static final String DFA5_acceptS =
        "\1\uffff\1\1\1\2\1\3\1\4\1\5\1\6\1\7\1\10\1\11\1\12\1\13";
    static final String DFA5_specialS =
        "\14\uffff}>";
    static final String[] DFA5_transitionS = {
            "\2\10\2\uffff\1\10\22\uffff\1\10\2\uffff\1\6\3\uffff\1\13\1"+
            "\2\1\4\2\uffff\1\3\1\uffff\1\11\1\uffff\12\11\1\5\1\1\37\uffff"+
            "\1\7",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
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
            return "1:1: Tokens : ( T__9 | T__10 | T__11 | T__12 | T__13 | T__14 | NESTED_ML_COMMENT | WS | DECIMAL_NUMBER | UNQUOTED_ALPHA_TEXT | QUOTED_TEXT );";
        }
    }
 

}