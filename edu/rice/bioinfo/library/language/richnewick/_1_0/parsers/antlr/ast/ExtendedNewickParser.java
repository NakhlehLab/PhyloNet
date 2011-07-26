// $ANTLR 3.3 Nov 30, 2010 12:45:30 C:\\Users\\Matt\\Desktop\\ExtendedNewick.g 2011-07-21 15:22:27

package edu.rice.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast;


import org.antlr.runtime.*;

public class ExtendedNewickParser extends Parser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "UNQUOTED_ALPHA_TEXT", "DECIMAL_NUMBER", "QUOTED_TEXT", "NESTED_ML_COMMENT", "WS", "';'", "'('", "','", "')'", "':'", "'#'"
    };
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


        public ExtendedNewickParser(TokenStream input) {
            this(input, new RecognizerSharedState());
        }
        public ExtendedNewickParser(TokenStream input, RecognizerSharedState state) {
            super(input, state);
             
        }
        

    public String[] getTokenNames() { return ExtendedNewickParser.tokenNames; }
    public String getGrammarFileName() { return "C:\\Users\\Matt\\Desktop\\ExtendedNewick.g"; }


    ParseStack stack = new ParseStackAction();

    ParseStack getParseStack()
    {
    	return stack;
    }

    protected void mismatch(IntStream input, int ttype, BitSet follow) throws RecognitionException
    {
    	throw new MismatchedTokenException(ttype, input);
    }

    public Object recoverFromMismatchedSet(IntStream input, RecognitionException e, BitSet follow) throws RecognitionException
    {
    	throw e;	
    }




    // $ANTLR start "network"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:43:1: network : (dl= descendant_list )? network_info ';' ;
    public final void network() throws RecognitionException {
        ExtendedNewickParser.descendant_list_return dl = null;


        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:43:9: ( (dl= descendant_list )? network_info ';' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:43:11: (dl= descendant_list )? network_info ';'
            {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:43:13: (dl= descendant_list )?
            int alt1=2;
            int LA1_0 = input.LA(1);

            if ( (LA1_0==10) ) {
                alt1=1;
            }
            switch (alt1) {
                case 1 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:43:13: dl= descendant_list
                    {
                    pushFollow(FOLLOW_descendant_list_in_network42);
                    dl=descendant_list();

                    state._fsp--;


                    }
                    break;

            }

            pushFollow(FOLLOW_network_info_in_network45);
            network_info();

            state._fsp--;

            match(input,9,FOLLOW_9_in_network47); 
             stack.pushNetwork((dl!=null?((Token)dl.start):null) != null); 

            }

        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "network"

    public static class descendant_list_return extends ParserRuleReturnScope {
    };

    // $ANTLR start "descendant_list"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:45:1: descendant_list : '(' subtree ( ',' subtree )* ')' ;
    public final ExtendedNewickParser.descendant_list_return descendant_list() throws RecognitionException {
        ExtendedNewickParser.descendant_list_return retval = new ExtendedNewickParser.descendant_list_return();
        retval.start = input.LT(1);

         int subTreeCount = 1; 
        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:45:48: ( '(' subtree ( ',' subtree )* ')' )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:45:50: '(' subtree ( ',' subtree )* ')'
            {
            match(input,10,FOLLOW_10_in_descendant_list63); 
            pushFollow(FOLLOW_subtree_in_descendant_list65);
            subtree();

            state._fsp--;

            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:45:62: ( ',' subtree )*
            loop2:
            do {
                int alt2=2;
                int LA2_0 = input.LA(1);

                if ( (LA2_0==11) ) {
                    alt2=1;
                }


                switch (alt2) {
            	case 1 :
            	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:45:63: ',' subtree
            	    {
            	    match(input,11,FOLLOW_11_in_descendant_list68); 
            	    pushFollow(FOLLOW_subtree_in_descendant_list70);
            	    subtree();

            	    state._fsp--;

            	    subTreeCount++;

            	    }
            	    break;

            	default :
            	    break loop2;
                }
            } while (true);

            match(input,12,FOLLOW_12_in_descendant_list77); 
             stack.pushDescendantList(subTreeCount); 

            }

            retval.stop = input.LT(-1);

        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return retval;
    }
    // $ANTLR end "descendant_list"


    // $ANTLR start "subtree"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:47:1: subtree : ( descendant_list network_info | network_info );
    public final void subtree() throws RecognitionException {
        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:47:9: ( descendant_list network_info | network_info )
            int alt3=2;
            int LA3_0 = input.LA(1);

            if ( (LA3_0==10) ) {
                alt3=1;
            }
            else if ( ((LA3_0>=UNQUOTED_ALPHA_TEXT && LA3_0<=QUOTED_TEXT)||LA3_0==9||(LA3_0>=11 && LA3_0<=14)) ) {
                alt3=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 3, 0, input);

                throw nvae;
            }
            switch (alt3) {
                case 1 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:47:11: descendant_list network_info
                    {
                    pushFollow(FOLLOW_descendant_list_in_subtree87);
                    descendant_list();

                    state._fsp--;

                    pushFollow(FOLLOW_network_info_in_subtree89);
                    network_info();

                    state._fsp--;

                     stack.pushSubtree(true); 

                    }
                    break;
                case 2 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:48:11: network_info
                    {
                    pushFollow(FOLLOW_network_info_in_subtree103);
                    network_info();

                    state._fsp--;

                     stack.pushSubtree(false); 

                    }
                    break;

            }
        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "subtree"


    // $ANTLR start "network_info"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:50:1: network_info : ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length bootstrap | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length bootstrap probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' bootstrap ( probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability );
    public final void network_info() throws RecognitionException {
        ExtendedNewickParser.node_label_return nl = null;

        ExtendedNewickParser.hybrid_node_qualifier_return hn = null;

        ExtendedNewickParser.branch_length_return bl = null;


        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:51:2: ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length bootstrap | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length bootstrap probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' bootstrap ( probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability )
            int alt16=5;
            alt16 = dfa16.predict(input);
            switch (alt16) {
                case 1 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:51:4: (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )?
                    {
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:51:6: (nl= node_label )?
                    int alt4=2;
                    int LA4_0 = input.LA(1);

                    if ( ((LA4_0>=UNQUOTED_ALPHA_TEXT && LA4_0<=QUOTED_TEXT)) ) {
                        alt4=1;
                    }
                    switch (alt4) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:51:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info132);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }

                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:51:21: (hn= hybrid_node_qualifier )?
                    int alt5=2;
                    int LA5_0 = input.LA(1);

                    if ( (LA5_0==14) ) {
                        alt5=1;
                    }
                    switch (alt5) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:51:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info137);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }

                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:51:47: (bl= branch_length )?
                    int alt6=2;
                    int LA6_0 = input.LA(1);

                    if ( (LA6_0==13) ) {
                        alt6=1;
                    }
                    switch (alt6) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:51:47: bl= branch_length
                            {
                            pushFollow(FOLLOW_branch_length_in_network_info142);
                            bl=branch_length();

                            state._fsp--;


                            }
                            break;

                    }

                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, (bl!=null?((Token)bl.start):null) != null, false, false); 

                    }
                    break;
                case 2 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:52:4: (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length bootstrap
                    {
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:52:6: (nl= node_label )?
                    int alt7=2;
                    int LA7_0 = input.LA(1);

                    if ( ((LA7_0>=UNQUOTED_ALPHA_TEXT && LA7_0<=QUOTED_TEXT)) ) {
                        alt7=1;
                    }
                    switch (alt7) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:52:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info171);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }

                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:52:21: (hn= hybrid_node_qualifier )?
                    int alt8=2;
                    int LA8_0 = input.LA(1);

                    if ( (LA8_0==14) ) {
                        alt8=1;
                    }
                    switch (alt8) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:52:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info176);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }

                    pushFollow(FOLLOW_branch_length_in_network_info182);
                    branch_length();

                    state._fsp--;

                    pushFollow(FOLLOW_bootstrap_in_network_info185);
                    bootstrap();

                    state._fsp--;

                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, true, true, false); 

                    }
                    break;
                case 3 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:53:4: (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length bootstrap probability
                    {
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:53:6: (nl= node_label )?
                    int alt9=2;
                    int LA9_0 = input.LA(1);

                    if ( ((LA9_0>=UNQUOTED_ALPHA_TEXT && LA9_0<=QUOTED_TEXT)) ) {
                        alt9=1;
                    }
                    switch (alt9) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:53:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info203);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }

                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:53:21: (hn= hybrid_node_qualifier )?
                    int alt10=2;
                    int LA10_0 = input.LA(1);

                    if ( (LA10_0==14) ) {
                        alt10=1;
                    }
                    switch (alt10) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:53:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info208);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }

                    pushFollow(FOLLOW_branch_length_in_network_info211);
                    branch_length();

                    state._fsp--;

                    pushFollow(FOLLOW_bootstrap_in_network_info214);
                    bootstrap();

                    state._fsp--;

                    pushFollow(FOLLOW_probability_in_network_info216);
                    probability();

                    state._fsp--;

                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, true, true, true); 

                    }
                    break;
                case 4 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:54:4: (nl= node_label )? (hn= hybrid_node_qualifier )? ':' bootstrap ( probability )?
                    {
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:54:6: (nl= node_label )?
                    int alt11=2;
                    int LA11_0 = input.LA(1);

                    if ( ((LA11_0>=UNQUOTED_ALPHA_TEXT && LA11_0<=QUOTED_TEXT)) ) {
                        alt11=1;
                    }
                    switch (alt11) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:54:6: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info225);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }

                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:54:21: (hn= hybrid_node_qualifier )?
                    int alt12=2;
                    int LA12_0 = input.LA(1);

                    if ( (LA12_0==14) ) {
                        alt12=1;
                    }
                    switch (alt12) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:54:21: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info230);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }

                    match(input,13,FOLLOW_13_in_network_info233); 
                    pushFollow(FOLLOW_bootstrap_in_network_info236);
                    bootstrap();

                    state._fsp--;

                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:54:60: ( probability )?
                    int alt13=2;
                    int LA13_0 = input.LA(1);

                    if ( (LA13_0==13) ) {
                        alt13=1;
                    }
                    switch (alt13) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:54:60: probability
                            {
                            pushFollow(FOLLOW_probability_in_network_info238);
                            probability();

                            state._fsp--;


                            }
                            break;

                    }

                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, false, true, true); 

                    }
                    break;
                case 5 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:55:5: (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability
                    {
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:55:7: (nl= node_label )?
                    int alt14=2;
                    int LA14_0 = input.LA(1);

                    if ( ((LA14_0>=UNQUOTED_ALPHA_TEXT && LA14_0<=QUOTED_TEXT)) ) {
                        alt14=1;
                    }
                    switch (alt14) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:55:7: nl= node_label
                            {
                            pushFollow(FOLLOW_node_label_in_network_info255);
                            nl=node_label();

                            state._fsp--;


                            }
                            break;

                    }

                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:55:22: (hn= hybrid_node_qualifier )?
                    int alt15=2;
                    int LA15_0 = input.LA(1);

                    if ( (LA15_0==14) ) {
                        alt15=1;
                    }
                    switch (alt15) {
                        case 1 :
                            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:55:22: hn= hybrid_node_qualifier
                            {
                            pushFollow(FOLLOW_hybrid_node_qualifier_in_network_info260);
                            hn=hybrid_node_qualifier();

                            state._fsp--;


                            }
                            break;

                    }

                    match(input,13,FOLLOW_13_in_network_info263); 
                    match(input,13,FOLLOW_13_in_network_info266); 
                    pushFollow(FOLLOW_probability_in_network_info268);
                    probability();

                    state._fsp--;

                     stack.pushNetworkInfo((nl!=null?((Token)nl.start):null) != null, (hn!=null?((Token)hn.start):null) != null, false, false, true); 

                    }
                    break;

            }
        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "network_info"

    public static class branch_length_return extends ParserRuleReturnScope {
    };

    // $ANTLR start "branch_length"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:57:1: branch_length : edge_label ;
    public final ExtendedNewickParser.branch_length_return branch_length() throws RecognitionException {
        ExtendedNewickParser.branch_length_return retval = new ExtendedNewickParser.branch_length_return();
        retval.start = input.LT(1);

        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:58:2: ( edge_label )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:58:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_branch_length289);
            edge_label();

            state._fsp--;

             stack.pushBranchLength(); 

            }

            retval.stop = input.LT(-1);

        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return retval;
    }
    // $ANTLR end "branch_length"


    // $ANTLR start "bootstrap"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:60:1: bootstrap : edge_label ;
    public final void bootstrap() throws RecognitionException {
        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:61:2: ( edge_label )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:61:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_bootstrap301);
            edge_label();

            state._fsp--;

             stack.pushBootstrap(); 

            }

        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "bootstrap"


    // $ANTLR start "probability"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:63:1: probability : edge_label ;
    public final void probability() throws RecognitionException {
        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:64:2: ( edge_label )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:64:4: edge_label
            {
            pushFollow(FOLLOW_edge_label_in_probability313);
            edge_label();

            state._fsp--;

             stack.pushProbability(); 

            }

        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "probability"

    public static class node_label_return extends ParserRuleReturnScope {
    };

    // $ANTLR start "node_label"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:66:1: node_label : text ;
    public final ExtendedNewickParser.node_label_return node_label() throws RecognitionException {
        ExtendedNewickParser.node_label_return retval = new ExtendedNewickParser.node_label_return();
        retval.start = input.LT(1);

        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:67:2: ( text )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:67:4: text
            {
            pushFollow(FOLLOW_text_in_node_label324);
            text();

            state._fsp--;

             stack.pushNodeLabel(); 

            }

            retval.stop = input.LT(-1);

        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return retval;
    }
    // $ANTLR end "node_label"

    public static class hybrid_node_qualifier_return extends ParserRuleReturnScope {
    };

    // $ANTLR start "hybrid_node_qualifier"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:69:1: hybrid_node_qualifier : '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER ;
    public final ExtendedNewickParser.hybrid_node_qualifier_return hybrid_node_qualifier() throws RecognitionException {
        ExtendedNewickParser.hybrid_node_qualifier_return retval = new ExtendedNewickParser.hybrid_node_qualifier_return();
        retval.start = input.LT(1);

        Token type=null;
        Token hybridNodeIndex=null;

        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:70:2: ( '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:70:4: '#' (type= UNQUOTED_ALPHA_TEXT )? hybridNodeIndex= DECIMAL_NUMBER
            {
            match(input,14,FOLLOW_14_in_hybrid_node_qualifier337); 
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:70:12: (type= UNQUOTED_ALPHA_TEXT )?
            int alt17=2;
            int LA17_0 = input.LA(1);

            if ( (LA17_0==UNQUOTED_ALPHA_TEXT) ) {
                alt17=1;
            }
            switch (alt17) {
                case 1 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:70:12: type= UNQUOTED_ALPHA_TEXT
                    {
                    type=(Token)match(input,UNQUOTED_ALPHA_TEXT,FOLLOW_UNQUOTED_ALPHA_TEXT_in_hybrid_node_qualifier341); 

                    }
                    break;

            }

            hybridNodeIndex=(Token)match(input,DECIMAL_NUMBER,FOLLOW_DECIMAL_NUMBER_in_hybrid_node_qualifier346); 
             stack.pushUnquotedText((hybridNodeIndex!=null?hybridNodeIndex.getText():null)); if(type!=null)stack.pushUnquotedText((type!=null?type.getText():null)); stack.pushHybridNodeQualifier(type, hybridNodeIndex); 

            }

            retval.stop = input.LT(-1);

        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return retval;
    }
    // $ANTLR end "hybrid_node_qualifier"


    // $ANTLR start "edge_label"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:73:1: edge_label : ':' text ;
    public final void edge_label() throws RecognitionException {
        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:74:2: ( ':' text )
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:74:4: ':' text
            {
            match(input,13,FOLLOW_13_in_edge_label360); 
            pushFollow(FOLLOW_text_in_edge_label362);
            text();

            state._fsp--;


            }

        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "edge_label"


    // $ANTLR start "text"
    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:77:1: text : (q= QUOTED_TEXT | (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+ );
    public final void text() throws RecognitionException {
        Token q=null;
        Token t=null;

         String str = ""; 
        try {
            // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:77:32: (q= QUOTED_TEXT | (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+ )
            int alt19=2;
            int LA19_0 = input.LA(1);

            if ( (LA19_0==QUOTED_TEXT) ) {
                alt19=1;
            }
            else if ( ((LA19_0>=UNQUOTED_ALPHA_TEXT && LA19_0<=DECIMAL_NUMBER)) ) {
                alt19=2;
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 19, 0, input);

                throw nvae;
            }
            switch (alt19) {
                case 1 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:77:34: q= QUOTED_TEXT
                    {
                    q=(Token)match(input,QUOTED_TEXT,FOLLOW_QUOTED_TEXT_in_text378); 
                     stack.pushQuotedText(q); 

                    }
                    break;
                case 2 :
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:78:7: (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+
                    {
                    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:78:7: (t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER ) )+
                    int cnt18=0;
                    loop18:
                    do {
                        int alt18=2;
                        int LA18_0 = input.LA(1);

                        if ( ((LA18_0>=UNQUOTED_ALPHA_TEXT && LA18_0<=DECIMAL_NUMBER)) ) {
                            alt18=1;
                        }


                        switch (alt18) {
                    	case 1 :
                    	    // C:\\Users\\Matt\\Desktop\\ExtendedNewick.g:78:8: t= ( UNQUOTED_ALPHA_TEXT | DECIMAL_NUMBER )
                    	    {
                    	    t=(Token)input.LT(1);
                    	    if ( (input.LA(1)>=UNQUOTED_ALPHA_TEXT && input.LA(1)<=DECIMAL_NUMBER) ) {
                    	        input.consume();
                    	        state.errorRecovery=false;
                    	    }
                    	    else {
                    	        MismatchedSetException mse = new MismatchedSetException(null,input);
                    	        throw mse;
                    	    }

                    	    str+=(t!=null?t.getText():null);

                    	    }
                    	    break;

                    	default :
                    	    if ( cnt18 >= 1 ) break loop18;
                                EarlyExitException eee =
                                    new EarlyExitException(18, input);
                                throw eee;
                        }
                        cnt18++;
                    } while (true);

                     stack.pushUnquotedText(str); 

                    }
                    break;

            }
        }

        catch(RecognitionException e)
        {
        	
        	throw e;
        }
        finally {
        }
        return ;
    }
    // $ANTLR end "text"

    // Delegated rules


    protected DFA16 dfa16 = new DFA16(this);
    static final String DFA16_eotS =
        "\22\uffff";
    static final String DFA16_eofS =
        "\22\uffff";
    static final String DFA16_minS =
        "\1\4\1\11\3\4\1\uffff\1\5\1\11\1\4\1\11\1\4\2\uffff\1\4\1\11\1"+
        "\4\2\uffff";
    static final String DFA16_maxS =
        "\3\16\1\5\1\15\1\uffff\1\5\4\15\2\uffff\1\6\2\15\2\uffff";
    static final String DFA16_acceptS =
        "\5\uffff\1\1\5\uffff\1\5\1\4\3\uffff\1\2\1\3";
    static final String DFA16_specialS =
        "\22\uffff}>";
    static final String[] DFA16_transitionS = {
            "\2\2\1\1\2\uffff\1\5\1\uffff\2\5\1\4\1\3",
            "\1\5\1\uffff\2\5\1\4\1\3",
            "\2\2\3\uffff\1\5\1\uffff\2\5\1\4\1\3",
            "\1\6\1\7",
            "\2\12\1\11\6\uffff\1\10",
            "",
            "\1\7",
            "\1\5\1\uffff\2\5\1\4",
            "\3\14\6\uffff\1\13",
            "\1\5\1\uffff\2\5\1\15",
            "\2\12\3\uffff\1\5\1\uffff\2\5\1\15",
            "",
            "",
            "\2\17\1\16",
            "\1\20\1\uffff\2\20\1\21",
            "\2\17\3\uffff\1\20\1\uffff\2\20\1\21",
            "",
            ""
    };

    static final short[] DFA16_eot = DFA.unpackEncodedString(DFA16_eotS);
    static final short[] DFA16_eof = DFA.unpackEncodedString(DFA16_eofS);
    static final char[] DFA16_min = DFA.unpackEncodedStringToUnsignedChars(DFA16_minS);
    static final char[] DFA16_max = DFA.unpackEncodedStringToUnsignedChars(DFA16_maxS);
    static final short[] DFA16_accept = DFA.unpackEncodedString(DFA16_acceptS);
    static final short[] DFA16_special = DFA.unpackEncodedString(DFA16_specialS);
    static final short[][] DFA16_transition;

    static {
        int numStates = DFA16_transitionS.length;
        DFA16_transition = new short[numStates][];
        for (int i=0; i<numStates; i++) {
            DFA16_transition[i] = DFA.unpackEncodedString(DFA16_transitionS[i]);
        }
    }

    class DFA16 extends DFA {

        public DFA16(BaseRecognizer recognizer) {
            this.recognizer = recognizer;
            this.decisionNumber = 16;
            this.eot = DFA16_eot;
            this.eof = DFA16_eof;
            this.min = DFA16_min;
            this.max = DFA16_max;
            this.accept = DFA16_accept;
            this.special = DFA16_special;
            this.transition = DFA16_transition;
        }
        public String getDescription() {
            return "50:1: network_info : ( (nl= node_label )? (hn= hybrid_node_qualifier )? (bl= branch_length )? | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length bootstrap | (nl= node_label )? (hn= hybrid_node_qualifier )? branch_length bootstrap probability | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' bootstrap ( probability )? | (nl= node_label )? (hn= hybrid_node_qualifier )? ':' ':' probability );";
        }
    }
 

    public static final BitSet FOLLOW_descendant_list_in_network42 = new BitSet(new long[]{0x0000000000006270L});
    public static final BitSet FOLLOW_network_info_in_network45 = new BitSet(new long[]{0x0000000000000200L});
    public static final BitSet FOLLOW_9_in_network47 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_10_in_descendant_list63 = new BitSet(new long[]{0x0000000000006470L});
    public static final BitSet FOLLOW_subtree_in_descendant_list65 = new BitSet(new long[]{0x0000000000001800L});
    public static final BitSet FOLLOW_11_in_descendant_list68 = new BitSet(new long[]{0x0000000000006470L});
    public static final BitSet FOLLOW_subtree_in_descendant_list70 = new BitSet(new long[]{0x0000000000001800L});
    public static final BitSet FOLLOW_12_in_descendant_list77 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_descendant_list_in_subtree87 = new BitSet(new long[]{0x0000000000006070L});
    public static final BitSet FOLLOW_network_info_in_subtree89 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_network_info_in_subtree103 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info132 = new BitSet(new long[]{0x0000000000006002L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info137 = new BitSet(new long[]{0x0000000000002002L});
    public static final BitSet FOLLOW_branch_length_in_network_info142 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info171 = new BitSet(new long[]{0x0000000000006000L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info176 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_branch_length_in_network_info182 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_bootstrap_in_network_info185 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info203 = new BitSet(new long[]{0x0000000000006000L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info208 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_branch_length_in_network_info211 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_bootstrap_in_network_info214 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_probability_in_network_info216 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info225 = new BitSet(new long[]{0x0000000000006000L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info230 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_13_in_network_info233 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_bootstrap_in_network_info236 = new BitSet(new long[]{0x0000000000002002L});
    public static final BitSet FOLLOW_probability_in_network_info238 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_node_label_in_network_info255 = new BitSet(new long[]{0x0000000000006000L});
    public static final BitSet FOLLOW_hybrid_node_qualifier_in_network_info260 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_13_in_network_info263 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_13_in_network_info266 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_probability_in_network_info268 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_branch_length289 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_bootstrap301 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_edge_label_in_probability313 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_text_in_node_label324 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_14_in_hybrid_node_qualifier337 = new BitSet(new long[]{0x0000000000000030L});
    public static final BitSet FOLLOW_UNQUOTED_ALPHA_TEXT_in_hybrid_node_qualifier341 = new BitSet(new long[]{0x0000000000000020L});
    public static final BitSet FOLLOW_DECIMAL_NUMBER_in_hybrid_node_qualifier346 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_13_in_edge_label360 = new BitSet(new long[]{0x0000000000000070L});
    public static final BitSet FOLLOW_text_in_edge_label362 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_QUOTED_TEXT_in_text378 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_text394 = new BitSet(new long[]{0x0000000000000032L});

}