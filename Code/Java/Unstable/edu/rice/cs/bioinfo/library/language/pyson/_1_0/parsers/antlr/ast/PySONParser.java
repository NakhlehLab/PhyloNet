// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2011-10-12 18:32:28

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;
import java.util.LinkedList;



import org.antlr.runtime.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONParser extends org.antlr.runtime.Parser {
    public static final String[] tokenNames = new String[] {
        "<invalid>", "<EOR>", "<DOWN>", "<UP>", "BEGIN", "DEFAULT_INDICATOR", "ELSE", "END", "ID", "ID_SET", "NESTED_ML_COMMENT", "NETWORK", "NETWORKS", "PHYLONET", "QUOTE", "ROOTAGE_QUALIFIER", "START", "TAXA_MAP", "TAXON_SET_LIST", "TRANSLATE", "TREE", "TREES", "UTREE", "WS", "'('", "')'", "','", "';'", "'='"
    };

    public static final int EOF=-1;
    public static final int T__24=24;
    public static final int T__25=25;
    public static final int T__26=26;
    public static final int T__27=27;
    public static final int T__28=28;
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
    public static final int TAXA_MAP=17;
    public static final int TAXON_SET_LIST=18;
    public static final int TRANSLATE=19;
    public static final int TREE=20;
    public static final int TREES=21;
    public static final int UTREE=22;
    public static final int WS=23;

    // delegates
    public Parser[] getDelegates() {
        return new Parser[] {};
    }

    // delegators


    public PySONParser(TokenStream input) {
        this(input, new RecognizerSharedState());
    }
    public PySONParser(TokenStream input, RecognizerSharedState state) {
        super(input, state);
    }

    public String[] getTokenNames() { return PySONParser.tokenNames; }
    public String getGrammarFileName() { return "D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g"; }



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

    @Override   
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





    // $ANTLR start "blocks"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:67:1: blocks : START ( block )* ;
    public final void blocks() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:67:9: ( START ( block )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:67:11: START ( block )*
            {
            match(input,START,FOLLOW_START_in_blocks51); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:67:17: ( block )*
            loop1:
            do {
                int alt1=2;
                int LA1_0 = input.LA(1);

                if ( (LA1_0==BEGIN) ) {
                    alt1=1;
                }


                switch (alt1) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:67:17: block
            	    {
            	    pushFollow(FOLLOW_block_in_blocks53);
            	    block();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop1;
                }
            } while (true);


             stack.pushBlocks(); 

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
    // $ANTLR end "blocks"



    // $ANTLR start "block"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:69:1: block : ( BEGIN networks_block_body end_semi | BEGIN trees_block_body end_semi | BEGIN phylonet_block_body end_semi | BEGIN skip_body end_semi );
    public final void block() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:69:7: ( BEGIN networks_block_body end_semi | BEGIN trees_block_body end_semi | BEGIN phylonet_block_body end_semi | BEGIN skip_body end_semi )
            int alt2=4;
            int LA2_0 = input.LA(1);

            if ( (LA2_0==BEGIN) ) {
                switch ( input.LA(2) ) {
                case NETWORKS:
                    {
                    alt2=1;
                    }
                    break;
                case TREES:
                    {
                    alt2=2;
                    }
                    break;
                case PHYLONET:
                    {
                    alt2=3;
                    }
                    break;
                case BEGIN:
                case DEFAULT_INDICATOR:
                case ELSE:
                case END:
                case ID:
                case ID_SET:
                case NESTED_ML_COMMENT:
                case NETWORK:
                case QUOTE:
                case ROOTAGE_QUALIFIER:
                case START:
                case TAXA_MAP:
                case TAXON_SET_LIST:
                case TRANSLATE:
                case TREE:
                case UTREE:
                case WS:
                case 24:
                case 25:
                case 26:
                case 27:
                case 28:
                    {
                    alt2=4;
                    }
                    break;
                default:
                    NoViableAltException nvae =
                        new NoViableAltException("", 2, 1, input);

                    throw nvae;

                }

            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 2, 0, input);

                throw nvae;

            }
            switch (alt2) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:69:9: BEGIN networks_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block64); 

                    pushFollow(FOLLOW_networks_block_body_in_block66);
                    networks_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block68);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:70:4: BEGIN trees_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block73); 

                    pushFollow(FOLLOW_trees_block_body_in_block75);
                    trees_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block80);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:71:4: BEGIN phylonet_block_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block85); 

                    pushFollow(FOLLOW_phylonet_block_body_in_block87);
                    phylonet_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block89);
                    end_semi();

                    state._fsp--;


                    }
                    break;
                case 4 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:72:4: BEGIN skip_body end_semi
                    {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block94); 

                    pushFollow(FOLLOW_skip_body_in_block96);
                    skip_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block108);
                    end_semi();

                    state._fsp--;


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
    // $ANTLR end "block"



    // $ANTLR start "networks_block_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:74:1: networks_block_body : ( NETWORKS ';' translation ( NETWORK rich_newick_assignment )* | NETWORKS ';' ( NETWORK rich_newick_assignment )* );
    public final void networks_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:75:2: ( NETWORKS ';' translation ( NETWORK rich_newick_assignment )* | NETWORKS ';' ( NETWORK rich_newick_assignment )* )
            int alt5=2;
            int LA5_0 = input.LA(1);

            if ( (LA5_0==NETWORKS) ) {
                int LA5_1 = input.LA(2);

                if ( (LA5_1==27) ) {
                    int LA5_2 = input.LA(3);

                    if ( (LA5_2==TRANSLATE) ) {
                        alt5=1;
                    }
                    else if ( (LA5_2==END||LA5_2==NETWORK) ) {
                        alt5=2;
                    }
                    else {
                        NoViableAltException nvae =
                            new NoViableAltException("", 5, 2, input);

                        throw nvae;

                    }
                }
                else {
                    NoViableAltException nvae =
                        new NoViableAltException("", 5, 1, input);

                    throw nvae;

                }
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 5, 0, input);

                throw nvae;

            }
            switch (alt5) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:76:3: NETWORKS ';' translation ( NETWORK rich_newick_assignment )*
                    {
                    match(input,NETWORKS,FOLLOW_NETWORKS_in_networks_block_body119); 

                    match(input,27,FOLLOW_27_in_networks_block_body121); 

                    pushFollow(FOLLOW_translation_in_networks_block_body123);
                    translation();

                    state._fsp--;


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:76:28: ( NETWORK rich_newick_assignment )*
                    loop3:
                    do {
                        int alt3=2;
                        int LA3_0 = input.LA(1);

                        if ( (LA3_0==NETWORK) ) {
                            alt3=1;
                        }


                        switch (alt3) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:76:29: NETWORK rich_newick_assignment
                    	    {
                    	    match(input,NETWORK,FOLLOW_NETWORK_in_networks_block_body126); 

                    	    pushFollow(FOLLOW_rich_newick_assignment_in_networks_block_body128);
                    	    rich_newick_assignment();

                    	    state._fsp--;


                    	     stack.pushNetworkAssignment(); 

                    	    }
                    	    break;

                    	default :
                    	    break loop3;
                        }
                    } while (true);


                     stack.pushNetworksBlockBody(true); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:77:4: NETWORKS ';' ( NETWORK rich_newick_assignment )*
                    {
                    match(input,NETWORKS,FOLLOW_NETWORKS_in_networks_block_body140); 

                    match(input,27,FOLLOW_27_in_networks_block_body142); 

                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:77:29: ( NETWORK rich_newick_assignment )*
                    loop4:
                    do {
                        int alt4=2;
                        int LA4_0 = input.LA(1);

                        if ( (LA4_0==NETWORK) ) {
                            alt4=1;
                        }


                        switch (alt4) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:77:30: NETWORK rich_newick_assignment
                    	    {
                    	    match(input,NETWORK,FOLLOW_NETWORK_in_networks_block_body157); 

                    	    pushFollow(FOLLOW_rich_newick_assignment_in_networks_block_body159);
                    	    rich_newick_assignment();

                    	    state._fsp--;


                    	     stack.pushNetworkAssignment(); 

                    	    }
                    	    break;

                    	default :
                    	    break loop4;
                        }
                    } while (true);


                     stack.pushNetworksBlockBody(false); 

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
    // $ANTLR end "networks_block_body"



    // $ANTLR start "trees_block_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:79:1: trees_block_body : ( TREES ';' translation ( tree_assigment )* | TREES ';' ( tree_assigment )* );
    public final void trees_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:80:2: ( TREES ';' translation ( tree_assigment )* | TREES ';' ( tree_assigment )* )
            int alt8=2;
            int LA8_0 = input.LA(1);

            if ( (LA8_0==TREES) ) {
                int LA8_1 = input.LA(2);

                if ( (LA8_1==27) ) {
                    int LA8_2 = input.LA(3);

                    if ( (LA8_2==TRANSLATE) ) {
                        alt8=1;
                    }
                    else if ( (LA8_2==END||LA8_2==TREE||LA8_2==UTREE) ) {
                        alt8=2;
                    }
                    else {
                        NoViableAltException nvae =
                            new NoViableAltException("", 8, 2, input);

                        throw nvae;

                    }
                }
                else {
                    NoViableAltException nvae =
                        new NoViableAltException("", 8, 1, input);

                    throw nvae;

                }
            }
            else {
                NoViableAltException nvae =
                    new NoViableAltException("", 8, 0, input);

                throw nvae;

            }
            switch (alt8) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:80:4: TREES ';' translation ( tree_assigment )*
                    {
                    match(input,TREES,FOLLOW_TREES_in_trees_block_body178); 

                    match(input,27,FOLLOW_27_in_trees_block_body180); 

                    pushFollow(FOLLOW_translation_in_trees_block_body182);
                    translation();

                    state._fsp--;


                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:80:27: ( tree_assigment )*
                    loop6:
                    do {
                        int alt6=2;
                        int LA6_0 = input.LA(1);

                        if ( (LA6_0==TREE||LA6_0==UTREE) ) {
                            alt6=1;
                        }


                        switch (alt6) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:80:27: tree_assigment
                    	    {
                    	    pushFollow(FOLLOW_tree_assigment_in_trees_block_body185);
                    	    tree_assigment();

                    	    state._fsp--;


                    	    }
                    	    break;

                    	default :
                    	    break loop6;
                        }
                    } while (true);


                     stack.pushTreesBlockBody(true); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:81:4: TREES ';' ( tree_assigment )*
                    {
                    match(input,TREES,FOLLOW_TREES_in_trees_block_body193); 

                    match(input,27,FOLLOW_27_in_trees_block_body195); 

                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:81:27: ( tree_assigment )*
                    loop7:
                    do {
                        int alt7=2;
                        int LA7_0 = input.LA(1);

                        if ( (LA7_0==TREE||LA7_0==UTREE) ) {
                            alt7=1;
                        }


                        switch (alt7) {
                    	case 1 :
                    	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:81:27: tree_assigment
                    	    {
                    	    pushFollow(FOLLOW_tree_assigment_in_trees_block_body210);
                    	    tree_assigment();

                    	    state._fsp--;


                    	    }
                    	    break;

                    	default :
                    	    break loop7;
                        }
                    } while (true);


                     stack.pushTreesBlockBody(false); 

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
    // $ANTLR end "trees_block_body"



    // $ANTLR start "phylonet_block_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:83:1: phylonet_block_body : PHYLONET ';' ( phylonet_command )* ;
    public final void phylonet_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:84:2: ( PHYLONET ';' ( phylonet_command )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:84:4: PHYLONET ';' ( phylonet_command )*
            {
            match(input,PHYLONET,FOLLOW_PHYLONET_in_phylonet_block_body222); 

            match(input,27,FOLLOW_27_in_phylonet_block_body224); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:84:17: ( phylonet_command )*
            loop9:
            do {
                int alt9=2;
                int LA9_0 = input.LA(1);

                if ( (LA9_0==BEGIN||(LA9_0 >= ID && LA9_0 <= ID_SET)||(LA9_0 >= NETWORK && LA9_0 <= QUOTE)||(LA9_0 >= TAXA_MAP && LA9_0 <= UTREE)||LA9_0==24||LA9_0==27) ) {
                    alt9=1;
                }


                switch (alt9) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:84:17: phylonet_command
            	    {
            	    pushFollow(FOLLOW_phylonet_command_in_phylonet_block_body226);
            	    phylonet_command();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop9;
                }
            } while (true);


             stack.pushPhylonetBlockBody(); 

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
    // $ANTLR end "phylonet_block_body"



    // $ANTLR start "phylonet_command"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:86:1: phylonet_command : ( phylonet_command_part )* ';' ;
    public final void phylonet_command() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:87:2: ( ( phylonet_command_part )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:87:5: ( phylonet_command_part )* ';'
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:87:5: ( phylonet_command_part )*
            loop10:
            do {
                int alt10=2;
                int LA10_0 = input.LA(1);

                if ( (LA10_0==BEGIN||(LA10_0 >= ID && LA10_0 <= ID_SET)||(LA10_0 >= NETWORK && LA10_0 <= QUOTE)||(LA10_0 >= TAXA_MAP && LA10_0 <= UTREE)||LA10_0==24) ) {
                    alt10=1;
                }


                switch (alt10) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:87:5: phylonet_command_part
            	    {
            	    pushFollow(FOLLOW_phylonet_command_part_in_phylonet_command240);
            	    phylonet_command_part();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop10;
                }
            } while (true);


            match(input,27,FOLLOW_27_in_phylonet_command243); 

             stack.pushPhylonetCommand(); 

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
    // $ANTLR end "phylonet_command"



    // $ANTLR start "phylonet_command_part"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:1: phylonet_command_part : ( identifier | ident_list |p= QUOTE |p= TAXON_SET_LIST |p= ID_SET |p= TAXA_MAP );
    public final void phylonet_command_part() throws RecognitionException {
        Token p=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:90:2: ( identifier | ident_list |p= QUOTE |p= TAXON_SET_LIST |p= ID_SET |p= TAXA_MAP )
            int alt11=6;
            switch ( input.LA(1) ) {
            case BEGIN:
            case ID:
            case NETWORK:
            case NETWORKS:
            case PHYLONET:
            case TRANSLATE:
            case TREE:
            case TREES:
            case UTREE:
                {
                alt11=1;
                }
                break;
            case 24:
                {
                alt11=2;
                }
                break;
            case QUOTE:
                {
                alt11=3;
                }
                break;
            case TAXON_SET_LIST:
                {
                alt11=4;
                }
                break;
            case ID_SET:
                {
                alt11=5;
                }
                break;
            case TAXA_MAP:
                {
                alt11=6;
                }
                break;
            default:
                NoViableAltException nvae =
                    new NoViableAltException("", 11, 0, input);

                throw nvae;

            }

            switch (alt11) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:90:4: identifier
                    {
                    pushFollow(FOLLOW_identifier_in_phylonet_command_part255);
                    identifier();

                    state._fsp--;


                     stack.pushPhylonetCommandPartIdent(); 

                    }
                    break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:4: ident_list
                    {
                    pushFollow(FOLLOW_ident_list_in_phylonet_command_part274);
                    ident_list();

                    state._fsp--;


                     stack.pushPhylonetCommandPartIdentList(); 

                    }
                    break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:92:4: p= QUOTE
                    {
                    p=(Token)match(input,QUOTE,FOLLOW_QUOTE_in_phylonet_command_part290); 

                     stack.pushPhylonetCommandPartQuote(p); 

                    }
                    break;
                case 4 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:93:10: p= TAXON_SET_LIST
                    {
                    p=(Token)match(input,TAXON_SET_LIST,FOLLOW_TAXON_SET_LIST_in_phylonet_command_part320); 

                     stack.pushPhylonetCommandPartSetList(p); 

                    }
                    break;
                case 5 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:94:4: p= ID_SET
                    {
                    p=(Token)match(input,ID_SET,FOLLOW_ID_SET_in_phylonet_command_part335); 

                     stack.pushPhylonetCommandPartIdSet(p); 

                    }
                    break;
                case 6 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:95:4: p= TAXA_MAP
                    {
                    p=(Token)match(input,TAXA_MAP,FOLLOW_TAXA_MAP_in_phylonet_command_part351); 

                     stack.pushPhylonetCommandPartTaxaMap(p); 

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
    // $ANTLR end "phylonet_command_part"



    // $ANTLR start "ident_list"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:97:1: ident_list : s= '(' identifier ( ',' identifier )* ')' ;
    public final void ident_list() throws RecognitionException {
        Token s=null;

         int numIdentsInList = 1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:2: (s= '(' identifier ( ',' identifier )* ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:4: s= '(' identifier ( ',' identifier )* ')'
            {
            s=(Token)match(input,24,FOLLOW_24_in_ident_list381); 

            pushFollow(FOLLOW_identifier_in_ident_list383);
            identifier();

            state._fsp--;


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:21: ( ',' identifier )*
            loop12:
            do {
                int alt12=2;
                int LA12_0 = input.LA(1);

                if ( (LA12_0==26) ) {
                    alt12=1;
                }


                switch (alt12) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:22: ',' identifier
            	    {
            	    match(input,26,FOLLOW_26_in_ident_list386); 

            	    pushFollow(FOLLOW_identifier_in_ident_list388);
            	    identifier();

            	    state._fsp--;


            	     numIdentsInList++; 

            	    }
            	    break;

            	default :
            	    break loop12;
                }
            } while (true);


            match(input,25,FOLLOW_25_in_ident_list395); 

             stack.pushIdentList(numIdentsInList, s); 

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
    // $ANTLR end "ident_list"



    // $ANTLR start "tree_assigment"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:100:1: tree_assigment : tr= ( TREE | UTREE ) rich_newick_assignment ;
    public final void tree_assigment() throws RecognitionException {
        Token tr=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:101:2: (tr= ( TREE | UTREE ) rich_newick_assignment )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:101:4: tr= ( TREE | UTREE ) rich_newick_assignment
            {
            tr=(Token)input.LT(1);

            if ( input.LA(1)==TREE||input.LA(1)==UTREE ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


            pushFollow(FOLLOW_rich_newick_assignment_in_tree_assigment417);
            rich_newick_assignment();

            state._fsp--;


             stack.pushTreeAssignment(tr); 

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
    // $ANTLR end "tree_assigment"



    // $ANTLR start "rich_newick_assignment"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:1: rich_newick_assignment : (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string ;
    public final void rich_newick_assignment() throws RecognitionException {
        Token d=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:104:2: ( (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:104:4: (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:104:5: (d= DEFAULT_INDICATOR )?
            int alt13=2;
            int LA13_0 = input.LA(1);

            if ( (LA13_0==DEFAULT_INDICATOR) ) {
                alt13=1;
            }
            switch (alt13) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:104:5: d= DEFAULT_INDICATOR
                    {
                    d=(Token)match(input,DEFAULT_INDICATOR,FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment431); 

                    }
                    break;

            }


            pushFollow(FOLLOW_identifier_in_rich_newick_assignment434);
            identifier();

            state._fsp--;


            match(input,28,FOLLOW_28_in_rich_newick_assignment436); 

            pushFollow(FOLLOW_rich_newick_string_in_rich_newick_assignment438);
            rich_newick_string();

            state._fsp--;


             stack.pushRichNewickAssignment(d==null); 

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
    // $ANTLR end "rich_newick_assignment"



    // $ANTLR start "translation"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:1: translation : TRANSLATE identifier WS identifier ( ',' identifier WS identifier )* ';' ;
    public final void translation() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:107:2: ( TRANSLATE identifier WS identifier ( ',' identifier WS identifier )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:107:4: TRANSLATE identifier WS identifier ( ',' identifier WS identifier )* ';'
            {
            match(input,TRANSLATE,FOLLOW_TRANSLATE_in_translation451); 

            pushFollow(FOLLOW_identifier_in_translation453);
            identifier();

            state._fsp--;


            match(input,WS,FOLLOW_WS_in_translation455); 

            pushFollow(FOLLOW_identifier_in_translation457);
            identifier();

            state._fsp--;


            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:107:39: ( ',' identifier WS identifier )*
            loop14:
            do {
                int alt14=2;
                int LA14_0 = input.LA(1);

                if ( (LA14_0==26) ) {
                    alt14=1;
                }


                switch (alt14) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:107:40: ',' identifier WS identifier
            	    {
            	    match(input,26,FOLLOW_26_in_translation460); 

            	    pushFollow(FOLLOW_identifier_in_translation462);
            	    identifier();

            	    state._fsp--;


            	    match(input,WS,FOLLOW_WS_in_translation464); 

            	    pushFollow(FOLLOW_identifier_in_translation466);
            	    identifier();

            	    state._fsp--;


            	    }
            	    break;

            	default :
            	    break loop14;
                }
            } while (true);


            match(input,27,FOLLOW_27_in_translation470); 

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
    // $ANTLR end "translation"



    // $ANTLR start "skip_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:109:1: skip_body : ~ ( TREES | NETWORKS | PHYLONET ) ';' (~ END )* ;
    public final void skip_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:110:2: (~ ( TREES | NETWORKS | PHYLONET ) ';' (~ END )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:110:4: ~ ( TREES | NETWORKS | PHYLONET ) ';' (~ END )*
            {
            if ( (input.LA(1) >= BEGIN && input.LA(1) <= NETWORK)||(input.LA(1) >= QUOTE && input.LA(1) <= TREE)||(input.LA(1) >= UTREE && input.LA(1) <= 28) ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


            match(input,27,FOLLOW_27_in_skip_body489); 

            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:110:35: (~ END )*
            loop15:
            do {
                int alt15=2;
                int LA15_0 = input.LA(1);

                if ( ((LA15_0 >= BEGIN && LA15_0 <= ELSE)||(LA15_0 >= ID && LA15_0 <= 28)) ) {
                    alt15=1;
                }


                switch (alt15) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
            	    {
            	    if ( (input.LA(1) >= BEGIN && input.LA(1) <= ELSE)||(input.LA(1) >= ID && input.LA(1) <= 28) ) {
            	        input.consume();
            	        state.errorRecovery=false;
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        throw mse;
            	    }


            	    }
            	    break;

            	default :
            	    break loop15;
                }
            } while (true);


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
    // $ANTLR end "skip_body"



    // $ANTLR start "rich_newick_string"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:113:1: rich_newick_string : (str=~ ( ';' ) )* ';' ;
    public final void rich_newick_string() throws RecognitionException {
        Token str=null;

         StringBuffer accum = new StringBuffer(); int line = -1; int col = -1; 
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:114:2: ( (str=~ ( ';' ) )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:114:5: (str=~ ( ';' ) )* ';'
            {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:114:5: (str=~ ( ';' ) )*
            loop16:
            do {
                int alt16=2;
                int LA16_0 = input.LA(1);

                if ( ((LA16_0 >= BEGIN && LA16_0 <= 26)||LA16_0==28) ) {
                    alt16=1;
                }


                switch (alt16) {
            	case 1 :
            	    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:114:6: str=~ ( ';' )
            	    {
            	    str=(Token)input.LT(1);

            	    if ( (input.LA(1) >= BEGIN && input.LA(1) <= 26)||input.LA(1)==28 ) {
            	        input.consume();
            	        state.errorRecovery=false;
            	    }
            	    else {
            	        MismatchedSetException mse = new MismatchedSetException(null,input);
            	        throw mse;
            	    }


            	    accum.append((str!=null?str.getText():null)); 
            	    	                      line = line == -1 ? (str!=null?str.getLine():0)    : line; 
            	    	                      col  = col  == -1 ? (str!=null?str.getCharPositionInLine():0)     : col ; 

            	    }
            	    break;

            	default :
            	    break loop16;
                }
            } while (true);


            match(input,27,FOLLOW_27_in_rich_newick_string535); 

             stack.pushRichNewickString(accum.toString() + ';', line, col); 

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
    // $ANTLR end "rich_newick_string"



    // $ANTLR start "identifier"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:119:1: identifier : s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | PHYLONET | ID ) ;
    public final void identifier() throws RecognitionException {
        Token s=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:2: (s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | PHYLONET | ID ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:4: s= ( TRANSLATE | TREE | UTREE | NETWORK | BEGIN | NETWORKS | TREES | PHYLONET | ID )
            {
            s=(Token)input.LT(1);

            if ( input.LA(1)==BEGIN||input.LA(1)==ID||(input.LA(1) >= NETWORK && input.LA(1) <= PHYLONET)||(input.LA(1) >= TRANSLATE && input.LA(1) <= UTREE) ) {
                input.consume();
                state.errorRecovery=false;
            }
            else {
                MismatchedSetException mse = new MismatchedSetException(null,input);
                throw mse;
            }


             stack.pushIdentifier(s); 

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
    // $ANTLR end "identifier"



    // $ANTLR start "end_semi"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:1: end_semi : END ';' ;
    public final void end_semi() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:9: ( END ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:11: END ';'
            {
            match(input,END,FOLLOW_END_in_end_semi594); 

            match(input,27,FOLLOW_27_in_end_semi596); 

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
    // $ANTLR end "end_semi"

    // Delegated rules


 

    public static final BitSet FOLLOW_START_in_blocks51 = new BitSet(new long[]{0x0000000000000012L});
    public static final BitSet FOLLOW_block_in_blocks53 = new BitSet(new long[]{0x0000000000000012L});
    public static final BitSet FOLLOW_BEGIN_in_block64 = new BitSet(new long[]{0x0000000000001000L});
    public static final BitSet FOLLOW_networks_block_body_in_block66 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_end_semi_in_block68 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block73 = new BitSet(new long[]{0x0000000000200000L});
    public static final BitSet FOLLOW_trees_block_body_in_block75 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_end_semi_in_block80 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block85 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_phylonet_block_body_in_block87 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_end_semi_in_block89 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block94 = new BitSet(new long[]{0x000000001FDFCFF0L});
    public static final BitSet FOLLOW_skip_body_in_block96 = new BitSet(new long[]{0x0000000000000080L});
    public static final BitSet FOLLOW_end_semi_in_block108 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_NETWORKS_in_networks_block_body119 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_networks_block_body121 = new BitSet(new long[]{0x0000000000080000L});
    public static final BitSet FOLLOW_translation_in_networks_block_body123 = new BitSet(new long[]{0x0000000000000802L});
    public static final BitSet FOLLOW_NETWORK_in_networks_block_body126 = new BitSet(new long[]{0x0000000000783930L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_networks_block_body128 = new BitSet(new long[]{0x0000000000000802L});
    public static final BitSet FOLLOW_NETWORKS_in_networks_block_body140 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_networks_block_body142 = new BitSet(new long[]{0x0000000000000802L});
    public static final BitSet FOLLOW_NETWORK_in_networks_block_body157 = new BitSet(new long[]{0x0000000000783930L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_networks_block_body159 = new BitSet(new long[]{0x0000000000000802L});
    public static final BitSet FOLLOW_TREES_in_trees_block_body178 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_trees_block_body180 = new BitSet(new long[]{0x0000000000080000L});
    public static final BitSet FOLLOW_translation_in_trees_block_body182 = new BitSet(new long[]{0x0000000000500002L});
    public static final BitSet FOLLOW_tree_assigment_in_trees_block_body185 = new BitSet(new long[]{0x0000000000500002L});
    public static final BitSet FOLLOW_TREES_in_trees_block_body193 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_trees_block_body195 = new BitSet(new long[]{0x0000000000500002L});
    public static final BitSet FOLLOW_tree_assigment_in_trees_block_body210 = new BitSet(new long[]{0x0000000000500002L});
    public static final BitSet FOLLOW_PHYLONET_in_phylonet_block_body222 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_phylonet_block_body224 = new BitSet(new long[]{0x00000000097E7B12L});
    public static final BitSet FOLLOW_phylonet_command_in_phylonet_block_body226 = new BitSet(new long[]{0x00000000097E7B12L});
    public static final BitSet FOLLOW_phylonet_command_part_in_phylonet_command240 = new BitSet(new long[]{0x00000000097E7B10L});
    public static final BitSet FOLLOW_27_in_phylonet_command243 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_phylonet_command_part255 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_ident_list_in_phylonet_command_part274 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_QUOTE_in_phylonet_command_part290 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_TAXON_SET_LIST_in_phylonet_command_part320 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_ID_SET_in_phylonet_command_part335 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_TAXA_MAP_in_phylonet_command_part351 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_24_in_ident_list381 = new BitSet(new long[]{0x0000000000783910L});
    public static final BitSet FOLLOW_identifier_in_ident_list383 = new BitSet(new long[]{0x0000000006000000L});
    public static final BitSet FOLLOW_26_in_ident_list386 = new BitSet(new long[]{0x0000000000783910L});
    public static final BitSet FOLLOW_identifier_in_ident_list388 = new BitSet(new long[]{0x0000000006000000L});
    public static final BitSet FOLLOW_25_in_ident_list395 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_tree_assigment409 = new BitSet(new long[]{0x0000000000783930L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_tree_assigment417 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment431 = new BitSet(new long[]{0x0000000000783910L});
    public static final BitSet FOLLOW_identifier_in_rich_newick_assignment434 = new BitSet(new long[]{0x0000000010000000L});
    public static final BitSet FOLLOW_28_in_rich_newick_assignment436 = new BitSet(new long[]{0x000000001FFFFFF0L});
    public static final BitSet FOLLOW_rich_newick_string_in_rich_newick_assignment438 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_TRANSLATE_in_translation451 = new BitSet(new long[]{0x0000000000783910L});
    public static final BitSet FOLLOW_identifier_in_translation453 = new BitSet(new long[]{0x0000000000800000L});
    public static final BitSet FOLLOW_WS_in_translation455 = new BitSet(new long[]{0x0000000000783910L});
    public static final BitSet FOLLOW_identifier_in_translation457 = new BitSet(new long[]{0x000000000C000000L});
    public static final BitSet FOLLOW_26_in_translation460 = new BitSet(new long[]{0x0000000000783910L});
    public static final BitSet FOLLOW_identifier_in_translation462 = new BitSet(new long[]{0x0000000000800000L});
    public static final BitSet FOLLOW_WS_in_translation464 = new BitSet(new long[]{0x0000000000783910L});
    public static final BitSet FOLLOW_identifier_in_translation466 = new BitSet(new long[]{0x000000000C000000L});
    public static final BitSet FOLLOW_27_in_translation470 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_skip_body480 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_skip_body489 = new BitSet(new long[]{0x000000001FFFFF72L});
    public static final BitSet FOLLOW_set_in_rich_newick_string514 = new BitSet(new long[]{0x000000001FFFFFF0L});
    public static final BitSet FOLLOW_27_in_rich_newick_string535 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_identifier549 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_END_in_end_semi594 = new BitSet(new long[]{0x0000000008000000L});
    public static final BitSet FOLLOW_27_in_end_semi596 = new BitSet(new long[]{0x0000000000000002L});

}