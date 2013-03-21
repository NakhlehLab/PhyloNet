// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g 2013-03-13 17:29:57

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.parsers.antlr.ast;
import java.util.LinkedList;


import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class RichNewick_1_1Parser extends Parser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "DECIMAL_NUMBER", "NESTED_ML_COMMENT", "QUOTED_TEXT", "ROOTAGE_QUALIFIER", "TREE_PROB", "UNQUOTED_ALPHA_TEXT", "WS", "'#'", "'('", "')'", "','", "':'", "';'"
    };

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

    // delegates
    public Parser[] getDelegates() {
        return new Parser[] {};
    }

    // delegators


    public RichNewick_1_1Parser(TokenStream input) {
        this(input, new RecognizerSharedState());
    }
    public RichNewick_1_1Parser(TokenStream input, RecognizerSharedState state) {
        super(input, state);
    }

    public String[] getTokenNames() { return RichNewick_1_1Parser.tokenNames; }
    public String getGrammarFileName() { return "D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g"; }



    public class ErrorWrapper
    {
    	public final String Message;
    	public final int Line,Col;
    	
    	ErrorWrapper(String message, int line, int col)
    	{
    		Message = message;
    		Line = line;
    		Col = col;
    	}
    }

    private List<ErrorWrapper> errors = new LinkedList<ErrorWrapper>();

       
    public void displayRecognitionError(String[] tokenNames, RecognitionException e) 
    {
            errors.add(new ErrorWrapper(getErrorMessage(e, tokenNames), e.line, e.c));
    }
        
    public List<ErrorWrapper> getErrors() 
    {
            return errors;
    }

    ParseStack stack = new ParseStackAction();

    ParseStack getParseStack()
    {
    	return stack;
    }





    // $ANTLR start "networks"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:62:1: networks : ( network )* ;
    public final void networks() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:62:9: ( ( network )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:62:11: ( network )*
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:62:11: ( network )*
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==DECIMAL_NUMBER||(LA1_0 >= QUOTED_TEXT && LA1_0 <= UNQUOTED_ALPHA_TEXT)||(LA1_0 >= 11 && LA1_0 <= 16)) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:62:11: network
            	    {
            	    pushFollow(FOLLOW_network_in_networks42);
            	    network();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop1;
                }
            } while (true);


             stack.pushNetworks(); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "networks"



    // $ANTLR start "network"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:1: network : (tp= TREE_PROB )? (rq= ROOTAGE_QUALIFIER )? (dl= descendant_list )? network_info ';' ;
    public final void network() throws RecognitionException {
        Token tp=null;
        Token rq=null;
        RichNewick_1_1Parser.descendant_list_return dl =null;


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:9: ( (tp= TREE_PROB )? (rq= ROOTAGE_QUALIFIER )? (dl= descendant_list )? network_info ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:11: (tp= TREE_PROB )? (rq= ROOTAGE_QUALIFIER )? (dl= descendant_list )? network_info ';'
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:13: (tp= TREE_PROB )?
            int alt2=2;
            int LA2_0 = input.LA(1);

            if ( (LA2_0==TREE_PROB) ) {
                alt2=1;
            }
            switch (alt2) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:13: tp= TREE_PROB
                    {
                    tp=(Token)match(input,TREE_PROB,FOLLOW_TREE_PROB_in_network55); 

                    }
                    break;

            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:27: (rq= ROOTAGE_QUALIFIER )?
            int alt3=2;
            int LA3_0 = input.LA(1);

            if ( (LA3_0==ROOTAGE_QUALIFIER) ) {
                alt3=1;
            }
            switch (alt3) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:27: rq= ROOTAGE_QUALIFIER
                    {
                    rq=(Token)match(input,ROOTAGE_QUALIFIER,FOLLOW_ROOTAGE_QUALIFIER_in_network60); 

                    }
                    break;

            }


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:49: (dl= descendant_list )?
            int alt4=2;
            int LA4_0 = input.LA(1);

            if ( (LA4_0==12) ) {
                alt4=1;
            }
            switch (alt4) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:64:49: dl= descendant_list
                    {
                    pushFollow(FOLLOW_descendant_list_in_network65);
                    dl=descendant_list();

                    state._fsp--;


                    }
                    break;

            }


            pushFollow(FOLLOW_network_info_in_network68);
            network_info();

            state._fsp--;


            match(input,16,FOLLOW_16_in_network70); 

             stack.pushNetwork(tp, rq, (dl!=null?((Token)dl.start):null) != null); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "network"


    public static class descendant_list_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "descendant_list"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:66:1: descendant_list : '(' subtree ( ',' subtree )* ')' ;
    public final RichNewick_1_1Parser.descendant_list_return descendant_list() throws RecognitionException {
        RichNewick_1_1Parser.descendant_list_return retval = new RichNewick_1_1Parser.descendant_list_return();
        retval.start = input.LT(1);


         int subTreeCount = 1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:66:48: ( '(' subtree ( ',' subtree )* ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:66:50: '(' subtree ( ',' subtree )* ')'
            {
            match(input,12,FOLLOW_12_in_descendant_list86); 

            pushFollow(FOLLOW_subtree_in_descendant_list88);
            subtree();

            state._fsp--;


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:66:62: ( ',' subtree )*
            loop5:
            do {
                int alt5=2;
                int LA5_0 = input.LA(1);

                if ( (LA5_0==14) ) {
                    alt5=1;
                }


                switch (alt5) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:66:63: ',' subtree
            	    {
            	    match(input,14,FOLLOW_14_in_descendant_list91); 

            	    pushFollow(FOLLOW_subtree_in_descendant_list93);
            	    subtree();

            	    state._fsp--;


            	    subTreeCount++;

            	    }
            	    break;

            	default :
            	    break loop5;
                }
            } while (true);


            match(input,13,FOLLOW_13_in_descendant_list100); 

             stack.pushDescendantList(subTreeCount); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "descendant_list"



    // $ANTLR start "subtree"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:68:1: subtree : ( descendant_list network_info | network_info );
    public final void subtree() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:68:9: ( descendant_list network_info | network_info )
            int alt6=2;
            int LA6_0 = input.LA(1);

            if ( (LA6_0==12) ) {
                alt6=1;
            }
            else if ( (LA6_0==DECIMAL_NUMBER||LA6_0==QUOTED_TEXT||LA6_0==UNQUOTED_ALPHA_TEXT||LA6_0==11||(LA6_0 >= 13 && LA6_0 <= 16)) ) {
                alt6=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 6, 0, input);

                throw nvae;

            }
            switch (alt6) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:68:11: descendant_list network_info
                    {
                    pushFollow(FOLLOW_descendant_list_in_subtree110);
                    descendant_list();

                    state._fsp--;


                    pushFollow(FOLLOW_network_info_in_subtree112);
                    network_info();

                    state._fsp--;


                     stack.pushSubtree(true); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:69:11: network_info
                    {
                    pushFollow(FOLLOW_network_info_in_subtree126);
                    network_info();

                    state._fsp--;


                     stack.pushSubtree(false); 

                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "subtree"



    // $ANTLR start "network_info"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:71:1: network_info : ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support probability | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length ':' probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' support (pb= probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability );
    public final void network_info() throws RecognitionException {
        RichNewick_1_1Parser.node_label_return nl =null;

        RichNewick_1_1Parser.hybrid_node_qualifier_return hn =null;

        RichNewick_1_1Parser.branch_length_return bl =null;

        RichNewick_1_1Parser.probability_return pb =null;


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:72:2: ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support probability | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length ':' probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' support (pb= probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability )
            int alt21=6;
            alt21 = dfa21.predict(input);
            switch (alt21) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:72:4: (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )?
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:72:6: (nl= node_label )?
                    int alt7=2;
                    int LA7_0 = input.LA(1);

                    if ( (LA7_0==DECIMAL_NUMBER||LA7_0==QUOTED_TEXT||LA7_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt7=1;
                    }
                    switch (alt7) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:72:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info155);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:72:21: (hn= hybrid_node_qualifier )?
                    int alt8=2;
                    int LA8_0 = input.LA(1);

                    if ( (LA8_0==11) ) {
                        alt8=1;
                    }
                    switch (alt8) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:72:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info160);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:72:47: (bl= branch_length )?
                    int alt9=2;
                    int LA9_0 = input.LA(1);

                    if ( (LA9_0==15) ) {
                        alt9=1;
                    }
                    switch (alt9) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:72:47: bl= branch_length
                            {
                            pushFollow(FOLLOW_branch_length_in_network_info165);
                            bl=branch_length();

                            state._fsp--;


                            }
                            break;

                    }


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, (bl!=null?((Token)bl.start):null) != null, false, false); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:73:4: (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:73:6: (nl= node_label )?
                    int alt10=2;
                    int LA10_0 = input.LA(1);

                    if ( (LA10_0==DECIMAL_NUMBER||LA10_0==QUOTED_TEXT||LA10_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt10=1;
                    }
                    switch (alt10) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:73:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info194);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:73:21: (hn= hybrid_node_qualifier )?
                    int alt11=2;
                    int LA11_0 = input.LA(1);

                    if ( (LA11_0==11) ) {
                        alt11=1;
                    }
                    switch (alt11) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:73:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info199);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    pushFollow(FOLLOW_branch_length_in_network_info205);
                    branch_length();

                    state._fsp--;


                    pushFollow(FOLLOW_support_in_network_info208);
                    support();

                    state._fsp--;


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, true, true, false); 

                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:74:4: (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support probability
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:74:6: (nl= node_label )?
                    int alt12=2;
                    int LA12_0 = input.LA(1);

                    if ( (LA12_0==DECIMAL_NUMBER||LA12_0==QUOTED_TEXT||LA12_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt12=1;
                    }
                    switch (alt12) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:74:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info228);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:74:21: (hn= hybrid_node_qualifier )?
                    int alt13=2;
                    int LA13_0 = input.LA(1);

                    if ( (LA13_0==11) ) {
                        alt13=1;
                    }
                    switch (alt13) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:74:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info233);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    pushFollow(FOLLOW_branch_length_in_network_info236);
                    branch_length();

                    state._fsp--;


                    pushFollow(FOLLOW_support_in_network_info239);
                    support();

                    state._fsp--;


                    pushFollow(FOLLOW_probability_in_network_info241);
                    probability();

                    state._fsp--;


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, true, true, true); 

                    }
                    break;
                case 4 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:75:4: (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length ':' probability
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:75:6: (nl= node_label )?
                    int alt14=2;
                    int LA14_0 = input.LA(1);

                    if ( (LA14_0==DECIMAL_NUMBER||LA14_0==QUOTED_TEXT||LA14_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt14=1;
                    }
                    switch (alt14) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:75:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info252);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:75:21: (hn= hybrid_node_qualifier )?
                    int alt15=2;
                    int LA15_0 = input.LA(1);

                    if ( (LA15_0==11) ) {
                        alt15=1;
                    }
                    switch (alt15) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:75:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info257);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    pushFollow(FOLLOW_branch_length_in_network_info260);
                    branch_length();

                    state._fsp--;


                    match(input,15,FOLLOW_15_in_network_info263); 

                    pushFollow(FOLLOW_probability_in_network_info265);
                    probability();

                    state._fsp--;


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, true, false, true); 

                    }
                    break;
                case 5 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:76:4: (nl= node_label )? (hn= hybrid_node_qualifier )? ':' support (pb= probability )?
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:76:6: (nl= node_label )?
                    int alt16=2;
                    int LA16_0 = input.LA(1);

                    if ( (LA16_0==DECIMAL_NUMBER||LA16_0==QUOTED_TEXT||LA16_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt16=1;
                    }
                    switch (alt16) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:76:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info280);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:76:21: (hn= hybrid_node_qualifier )?
                    int alt17=2;
                    int LA17_0 = input.LA(1);

                    if ( (LA17_0==11) ) {
                        alt17=1;
                    }
                    switch (alt17) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:76:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info285);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    match(input,15,FOLLOW_15_in_network_info288); 

                    pushFollow(FOLLOW_support_in_network_info291);
                    support();

                    state._fsp--;


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:76:60: (pb= probability )?
                    int alt18=2;
                    int LA18_0 = input.LA(1);

                    if ( (LA18_0==15) ) {
                        alt18=1;
                    }
                    switch (alt18) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:76:60: pb= probability
                            {
                            pushFollow(FOLLOW_probability_in_network_info295);
                            pb=probability();

                            state._fsp--;


                            }
                            break;

                    }


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, false, true, (pb!=null?((Token)pb.start):null) != null); 

                    }
                    break;
                case 6 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:77:5: (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:77:7: (nl= node_label )?
                    int alt19=2;
                    int LA19_0 = input.LA(1);

                    if ( (LA19_0==DECIMAL_NUMBER||LA19_0==QUOTED_TEXT||LA19_0==UNQUOTED_ALPHA_TEXT) ) {
                        alt19=1;
                    }
                    switch (alt19) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:77:7: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info312);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:77:22: (hn= hybrid_node_qualifier )?
                    int alt20=2;
                    int LA20_0 = input.LA(1);

                    if ( (LA20_0==11) ) {
                        alt20=1;
                    }
                    switch (alt20) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:77:22: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info317);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }


                    match(input,15,FOLLOW_15_in_network_info320); 

                    match(input,15,FOLLOW_15_in_network_info323); 

                    pushFollow(FOLLOW_probability_in_network_info325);
                    probability();

                    state._fsp--;


                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, false, false, true); 

                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "network_info"


    public static class branch_length_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "branch_length"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:79:1: branch_length : edge_label ;
    public final RichNewick_1_1Parser.branch_length_return branch_length() throws RecognitionException {
        RichNewick_1_1Parser.branch_length_return retval = new RichNewick_1_1Parser.branch_length_return();
        retval.start = input.LT(1);


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:80:2: ( edge_label )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:80:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_branch_length346);
            edge_label();

            state._fsp--;


             stack.pushBranchLength(); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "branch_length"



    // $ANTLR start "support"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:82:1: support : edge_label ;
    public final void support() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:83:2: ( edge_label )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:83:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_support358);
            edge_label();

            state._fsp--;


             stack.pushSupport(); 

            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "support"


    public static class probability_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "probability"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:85:1: probability : edge_label ;
    public final RichNewick_1_1Parser.probability_return probability() throws RecognitionException {
        RichNewick_1_1Parser.probability_return retval = new RichNewick_1_1Parser.probability_return();
        retval.start = input.LT(1);


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:86:2: ( edge_label )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:86:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_probability370);
            edge_label();

            state._fsp--;


             stack.pushProbability(); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "probability"


    public static class node_label_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "node_label"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:88:1: node_label : text ;
    public final RichNewick_1_1Parser.node_label_return node_label() throws RecognitionException {
        RichNewick_1_1Parser.node_label_return retval = new RichNewick_1_1Parser.node_label_return();
        retval.start = input.LT(1);


        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:89:2: ( text )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:89:4: text
            {
            pushFollow(FOLLOW_text_in_node_label381);
            text();

            state._fsp--;


             stack.pushNodeLabel(); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "node_label"


    public static class hybrid_node_qualifier_return extends ParserRuleReturnScope {
    };


    // $ANTLR start "hybrid_node_qualifier"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:91:1: hybrid_node_qualifier : '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER ;
    public final RichNewick_1_1Parser.hybrid_node_qualifier_return hybrid_node_qualifier() throws RecognitionException {
        RichNewick_1_1Parser.hybrid_node_qualifier_return retval = new RichNewick_1_1Parser.hybrid_node_qualifier_return();
        retval.start = input.LT(1);


        Token type=null;
        Token hybridNodeIndex=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:92:2: ( '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:92:4: '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER
            {
            match(input,11,FOLLOW_11_in_hybrid_node_qualifier394); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:92:12: (type= UNQUOTED_ALPHA_TEXT )?
            int alt22=2;
            int LA22_0 = input.LA(1);

            if ( (LA22_0==UNQUOTED_ALPHA_TEXT) ) {
                alt22=1;
            }
            switch (alt22) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:92:12: type= UNQUOTED_ALPHA_TEXT
                    {
                    type=(Token)match(input,UNQUOTED_ALPHA_TEXT,FOLLOW_UNQUOTED_ALPHA_TEXT_in_hybrid_node_qualifier398); 

                    }
                    break;

            }


            hybridNodeIndex=(Token)match(input,DECIMAL_NUMBER,FOLLOW_DECIMAL_NUMBER_in_hybrid_node_qualifier403); 

             stack.pushUnquotedText((hybridNodeIndex!=null?hybridNodeIndex.getText():null), (hybridNodeIndex!=null?hybridNodeIndex.getLine():0), (hybridNodeIndex!=null?hybridNodeIndex.getCharPositionInLine():0)); 
            		  if(type!=null)stack.pushUnquotedText((type!=null?type.getText():null), (type!=null?type.getLine():0), (type!=null?type.getCharPositionInLine():0)); stack.pushHybridNodeQualifier(type, hybridNodeIndex); 

            }

            retval.stop = input.LT(-1);


        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return retval;
    }
    // $ANTLR end "hybrid_node_qualifier"



    // $ANTLR start "edge_label"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:96:1: edge_label : ':' text ;
    public final void edge_label() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:97:2: ( ':' text )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:97:4: ':' text
            {
            match(input,15,FOLLOW_15_in_edge_label417); 

            pushFollow(FOLLOW_text_in_edge_label419);
            text();

            state._fsp--;


            }

        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "edge_label"



    // $ANTLR start "text"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:100:1: text : (q= QUOTED_TEXT | (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+ );
    public final void text() throws RecognitionException {
        Token q=null;
        Token t=null;

         String str = ""; int ln = -1; int cn = -1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:100:58: (q= QUOTED_TEXT | (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+ )
            int alt24=2;
            int LA24_0 = input.LA(1);

            if ( (LA24_0==QUOTED_TEXT) ) {
                alt24=1;
            }
            else if ( (LA24_0==DECIMAL_NUMBER||LA24_0==UNQUOTED_ALPHA_TEXT) ) {
                alt24=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 24, 0, input);

                throw nvae;

            }
            switch (alt24) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:101:14: q= QUOTED_TEXT
                    {
                    q=(Token)match(input,QUOTED_TEXT,FOLLOW_QUOTED_TEXT_in_text449); 

                     stack.pushQuotedText(q); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:102:7: (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+
                    {
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:102:7: (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+
                    int cnt23=0;
                    loop23:
                    do {
                        int alt23=2;
                        int LA23_0 = input.LA(1);

                        if ( (LA23_0==DECIMAL_NUMBER||LA23_0==UNQUOTED_ALPHA_TEXT) ) {
                            alt23=1;
                        }


                        switch (alt23) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\RichNewick\\RichNewick_1_1.g:102:8: t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER )
                    	    {
                    	    t=(Token)input.LT(1);

                    	    if ( input.LA(1)==DECIMAL_NUMBER||input.LA(1)==UNQUOTED_ALPHA_TEXT ) {
                    	        input.consume();
                    	        state.errorRecovery=false;
                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        throw mse;
                    	    }


                    	    str+=(t!=null?t.getText():null); ln = ln == -1 ? t.getLine() : ln; cn = cn == -1 ? t.getCharPositionInLine() : cn;

                    	    }
                    	    break;

                    	default :
                    	    if ( cnt23 >= 1 ) break loop23;
                                EarlyExitException eee =
                                    new EarlyExitException(23, input);
                                throw eee;
                        }
                        cnt23++;
                    } while (true);


                     
                    				  

                            				str = str.trim();

                            				if(str.contains(" ") || str.contains("\n") || str.contains("\r") || str.contains("\t") )
                            				{
                                					errors.add(new ErrorWrapper("Invalid whitespace in node or edge label '" + str + "'.", ln, cn));
                            				}
                    				  
                    				  	stack.pushUnquotedText(str, ln, cn); 
                    				  

                    }
                    break;

            }
        }
        catch (RecognitionException re) {
            reportError(re);
            recover(input,re);
        }

        finally {
        	// do for sure before leaving
        }
        return ;
    }
    // $ANTLR end "text"

    // Delegated rules


    protected DFA21 dfa21 = new DFA21(this);
    static final String DFA21_eotS =
        "\23\uffff";
    static final String DFA21_eofS =
        "\23\uffff";
    static final String DFA21_minS =
        "\1\4\1\13\3\4\1\uffff\1\4\1\15\1\4\1\15\1\4\2\uffff\1\4\1\15\1\4"+
        "\3\uffff";
    static final String DFA21_maxS =
        "\3\20\1\11\1\17\1\uffff\1\4\1\20\1\17\2\20\2\uffff\1\17\2\20\3\uffff";
    static final String DFA21_acceptS =
        "\5\uffff\1\1\5\uffff\1\6\1\5\3\uffff\1\4\1\2\1\3";
    static final String DFA21_specialS =
        "\23\uffff}>";
    static final String[] DFA21_transitionS = {
            "\1\2\1\uffff\1\1\2\uffff\1\2\1\uffff\1\3\1\uffff\2\5\1\4\1\5",
            "\1\3\1\uffff\2\5\1\4\1\5",
            "\1\2\4\uffff\1\2\1\uffff\1\3\1\uffff\2\5\1\4\1\5",
            "\1\7\4\uffff\1\6",
            "\1\12\1\uffff\1\11\2\uffff\1\12\5\uffff\1\10",
            "",
            "\1\7",
            "\2\5\1\4\1\5",
            "\1\14\1\uffff\1\14\2\uffff\1\14\5\uffff\1\13",
            "\2\5\1\15\1\5",
            "\1\12\4\uffff\1\12\3\uffff\2\5\1\15\1\5",
            "",
            "",
            "\1\17\1\uffff\1\16\2\uffff\1\17\5\uffff\1\20",
            "\2\21\1\22\1\21",
            "\1\17\4\uffff\1\17\3\uffff\2\21\1\22\1\21",
            "",
            "",
            ""
    };

    static final short[] DFA21_eot = DFA.unpackEncodedString(DFA21_eotS);
    static final short[] DFA21_eof = DFA.unpackEncodedString(DFA21_eofS);
    static final char[] DFA21_min = DFA.unpackEncodedStringToUnsignedChars(DFA21_minS);
    static final char[] DFA21_max = DFA.unpackEncodedStringToUnsignedChars(DFA21_maxS);
    static final short[] DFA21_accept = DFA.unpackEncodedString(DFA21_acceptS);
    static final short[] DFA21_special = DFA.unpackEncodedString(DFA21_specialS);
    static final short[][] DFA21_transition;

    static {
        int numStates = DFA21_transitionS.length;
        DFA21_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA21_transition[i] = DFA.unpackEncodedString(DFA21_transitionS[i]);
        }
    }

    class DFA21 extends DFA {

        public DFA21(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 21;
            this.eot = DFA21_eot;
            this.eof = DFA21_eof;
            this.min = DFA21_min;
            this.max = DFA21_max;
            this.accept = DFA21_accept;
            this.special = DFA21_special;
            this.transition = DFA21_transition;
        }
        public String getDescription() {
            return "71:1: network_info : ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length support probability | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length ':' probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' support (pb= probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability );";
        }
    }
 

    public static final BitSet FOLLOW_network_in_networks42 = new BitSet(new long[]{0x0000000000019BD2L});
    public static final BitSet FOLLOW_TREE_PROB_in_network55 = new BitSet(new long[]{0x0000000000019AD0L});
    public static final BitSet FOLLOW_ROOTAGE_QUALIFIER_in_network60 = new BitSet(new long[]{0x0000000000019A50L});
    public static final BitSet FOLLOW_descendant_list_in_network65 = new BitSet(new long[]{0x0000000000018A50L});
    public static final BitSet FOLLOW_network_info_in_network68 = new BitSet(new long[]{0x0000000000010000L});
    public static final BitSet FOLLOW_16_in_network70 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_12_in_descendant_list86 = new BitSet(new long[]{0x0000000000009A50L});
    public static final BitSet FOLLOW_subtree_in_descendant_list88 = new BitSet(new long[]{0x0000000000006000L});
    public static final BitSet FOLLOW_14_in_descendant_list91 = new BitSet(new long[]{0x0000000000009A50L});
    public static final BitSet FOLLOW_subtree_in_descendant_list93 = new BitSet(new long[]{0x0000000000006000L});
    public static final BitSet FOLLOW_13_in_descendant_list100 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_descendant_list_in_subtree110 = new BitSet(new long[]{0x0000000000008A50L});
    public static final BitSet FOLLOW_network_info_in_subtree112 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_network_info_in_subtree126 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info155 = new BitSet(new long[]{0x0000000000008802L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info160 = new BitSet(new long[]{0x0000000000008002L});
    public static final BitSet FOLLOW_branch_length_in_network_info165 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info194 = new BitSet(new long[]{0x0000000000008800L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info199 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_branch_length_in_network_info205 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_support_in_network_info208 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info228 = new BitSet(new long[]{0x0000000000008800L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info233 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_branch_length_in_network_info236 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_support_in_network_info239 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_probability_in_network_info241 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info252 = new BitSet(new long[]{0x0000000000008800L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info257 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_branch_length_in_network_info260 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_15_in_network_info263 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_probability_in_network_info265 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info280 = new BitSet(new long[]{0x0000000000008800L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info285 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_15_in_network_info288 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_support_in_network_info291 = new BitSet(new long[]{0x0000000000008002L});
    public static final BitSet FOLLOW_probability_in_network_info295 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info312 = new BitSet(new long[]{0x0000000000008800L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info317 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_15_in_network_info320 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_15_in_network_info323 = new BitSet(new long[]{0x0000000000008000L});
    public static final BitSet FOLLOW_probability_in_network_info325 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_branch_length346 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_support358 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_probability370 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_text_in_node_label381 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_11_in_hybrid_node_qualifier394 = new BitSet(new long[]{0x0000000000000210L});
    public static final BitSet FOLLOW_UNQUOTED_ALPHA_TEXT_in_hybrid_node_qualifier398 = new BitSet(new long[]{0x0000000000000010L});
    public static final BitSet FOLLOW_DECIMAL_NUMBER_in_hybrid_node_qualifier403 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_15_in_edge_label417 = new BitSet(new long[]{0x0000000000000250L});
    public static final BitSet FOLLOW_text_in_edge_label419 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_QUOTED_TEXT_in_text449 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_text465 = new BitSet(new long[]{0x0000000000000212L});

}