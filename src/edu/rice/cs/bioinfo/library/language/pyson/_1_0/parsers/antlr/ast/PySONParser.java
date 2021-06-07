// $ANTLR 3.4 D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g 2013-03-05 15:54:00

package edu.rice.cs.bioinfo.library.language.pyson._1_0.parsers.antlr.ast;
import org.antlr.runtime.*;

import java.util.LinkedList;
import java.util.List;

@SuppressWarnings({"all", "warnings", "unchecked"})
public class PySONParser extends org.antlr.runtime.Parser {
    public static final String[] tokenNames = new String[] {
            "<invalid>", "<EOR>", "<DOWN>", "<UP>", "BEGIN", "DATA", "DATATYPE", "DEFAULT_INDICATOR", "DIMENSIONS", "ELSE", "END", "FORMAT", "GAP", "ID", "ID_SET", "MATRIX", "MISSING", "NCHAR", "NESTED_ML_COMMENT", "NETWORK", "NETWORKS", "NTAX", "PHYLONET", "QUOTE", "RN_LS_NONCOMMENT", "START", "SYMBOLS", "TAXON_SET_LIST", "TRANSLATE", "TREE", "TREES", "UTREE", "WS", "'('", "')'", "','", "':'", "';'", "'<'", "'='", "'>'"
    };

    public static final int EOF=-1;
    public static final int T__33=33;
    public static final int T__34=34;
    public static final int T__35=35;
    public static final int T__36=36;
    public static final int T__37=37;
    public static final int T__38=38;
    public static final int T__39=39;
    public static final int T__40=40;
    public static final int BEGIN=4;
    public static final int DATA=5;
    public static final int DATATYPE=6;
    public static final int DEFAULT_INDICATOR=7;
    public static final int DIMENSIONS=8;
    public static final int ELSE=9;
    public static final int END=10;
    public static final int FORMAT=11;
    public static final int GAP=12;
    public static final int ID=13;
    public static final int ID_SET=14;
    public static final int MATRIX=15;
    public static final int MISSING=16;
    public static final int NCHAR=17;
    public static final int NESTED_ML_COMMENT=18;
    public static final int NETWORK=19;
    public static final int NETWORKS=20;
    public static final int NTAX=21;
    public static final int PHYLONET=22;
    public static final int QUOTE=23;
    public static final int RN_LS_NONCOMMENT=24;
    public static final int START=25;
    public static final int SYMBOLS=26;
    public static final int TAXON_SET_LIST=27;
    public static final int TRANSLATE=28;
    public static final int TREE=29;
    public static final int TREES=30;
    public static final int UTREE=31;
    public static final int WS=32;

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:1: blocks : START ( block )* EOF ;
    public final void blocks() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:9: ( START ( block )* EOF )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:11: START ( block )* EOF
            {
                match(input,START,FOLLOW_START_in_blocks52);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:17: ( block )*
                loop1:
                do {
                    int alt1=2;
                    int LA1_0 = input.LA(1);

                    if ( (LA1_0==BEGIN) ) {
                        alt1=1;
                    }


                    switch (alt1) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:89:17: block
                        {
                            pushFollow(FOLLOW_block_in_blocks54);
                            block();

                            state._fsp--;


                        }
                        break;

                        default :
                            break loop1;
                    }
                } while (true);


                stack.pushBlocks();

                match(input,EOF,FOLLOW_EOF_in_blocks59);

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:1: block : ( BEGIN networks_block_body end_semi | BEGIN trees_block_body end_semi | BEGIN phylonet_block_body end_semi | BEGIN data_block_body end_semi | BEGIN skip_body end_semi );
    public final void block() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:7: ( BEGIN networks_block_body end_semi | BEGIN trees_block_body end_semi | BEGIN phylonet_block_body end_semi | BEGIN data_block_body end_semi | BEGIN skip_body end_semi )
            int alt2=5;
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
                    case DATA:
                    {
                        alt2=4;
                    }
                    break;
                    case BEGIN:
                    case DATATYPE:
                    case DEFAULT_INDICATOR:
                    case DIMENSIONS:
                    case ELSE:
                    case END:
                    case FORMAT:
                    case GAP:
                    case ID:
                    case ID_SET:
                    case MATRIX:
                    case MISSING:
                    case NCHAR:
                    case NESTED_ML_COMMENT:
                    case NETWORK:
                    case NTAX:
                    case QUOTE:
                    case RN_LS_NONCOMMENT:
                    case START:
                    case SYMBOLS:
                    case TAXON_SET_LIST:
                    case TRANSLATE:
                    case TREE:
                    case UTREE:
                    case WS:
                    case 33:
                    case 34:
                    case 35:
                    case 36:
                    case 37:
                    case 38:
                    case 39:
                    case 40:
                    {
                        alt2=5;
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
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:91:9: BEGIN networks_block_body end_semi
                {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block67);

                    pushFollow(FOLLOW_networks_block_body_in_block69);
                    networks_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block71);
                    end_semi();

                    state._fsp--;


                }
                break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:92:4: BEGIN trees_block_body end_semi
                {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block76);

                    pushFollow(FOLLOW_trees_block_body_in_block78);
                    trees_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block83);
                    end_semi();

                    state._fsp--;


                }
                break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:93:4: BEGIN phylonet_block_body end_semi
                {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block88);

                    pushFollow(FOLLOW_phylonet_block_body_in_block90);
                    phylonet_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block92);
                    end_semi();

                    state._fsp--;


                }
                break;
                case 4 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:94:4: BEGIN data_block_body end_semi
                {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block97);

                    pushFollow(FOLLOW_data_block_body_in_block99);
                    data_block_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block103);
                    end_semi();

                    state._fsp--;


                }
                break;
                case 5 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:95:4: BEGIN skip_body end_semi
                {
                    match(input,BEGIN,FOLLOW_BEGIN_in_block108);

                    pushFollow(FOLLOW_skip_body_in_block110);
                    skip_body();

                    state._fsp--;


                    pushFollow(FOLLOW_end_semi_in_block122);
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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:97:1: networks_block_body : NETWORKS ';' ( NETWORK rich_newick_assignment )* ;
    public final void networks_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:98:2: ( NETWORKS ';' ( NETWORK rich_newick_assignment )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:99:3: NETWORKS ';' ( NETWORK rich_newick_assignment )*
            {
                match(input,NETWORKS,FOLLOW_NETWORKS_in_networks_block_body133);

                match(input,37,FOLLOW_37_in_networks_block_body135);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:99:28: ( NETWORK rich_newick_assignment )*
                loop3:
                do {
                    int alt3=2;
                    int LA3_0 = input.LA(1);

                    if ( (LA3_0==NETWORK) ) {
                        alt3=1;
                    }


                    switch (alt3) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:99:29: NETWORK rich_newick_assignment
                        {
                            match(input,NETWORK,FOLLOW_NETWORK_in_networks_block_body150);

                            pushFollow(FOLLOW_rich_newick_assignment_in_networks_block_body152);
                            rich_newick_assignment();

                            state._fsp--;


                            stack.pushNetworkAssignment();

                        }
                        break;

                        default :
                            break loop3;
                    }
                } while (true);


                stack.pushNetworksBlockBody(false);

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:101:1: trees_block_body : TREES ';' ( tree_assigment )* ;
    public final void trees_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:2: ( TREES ';' ( tree_assigment )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:4: TREES ';' ( tree_assigment )*
            {
                match(input,TREES,FOLLOW_TREES_in_trees_block_body172);

                match(input,37,FOLLOW_37_in_trees_block_body174);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:27: ( tree_assigment )*
                loop4:
                do {
                    int alt4=2;
                    int LA4_0 = input.LA(1);

                    if ( (LA4_0==TREE||LA4_0==UTREE) ) {
                        alt4=1;
                    }


                    switch (alt4) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:103:27: tree_assigment
                        {
                            pushFollow(FOLLOW_tree_assigment_in_trees_block_body189);
                            tree_assigment();

                            state._fsp--;


                        }
                        break;

                        default :
                            break loop4;
                    }
                } while (true);


                stack.pushTreesBlockBody(false);

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:105:1: phylonet_block_body : PHYLONET ';' ( phylonet_command )* ;
    public final void phylonet_block_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:2: ( PHYLONET ';' ( phylonet_command )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:4: PHYLONET ';' ( phylonet_command )*
            {
                match(input,PHYLONET,FOLLOW_PHYLONET_in_phylonet_block_body201);

                match(input,37,FOLLOW_37_in_phylonet_block_body203);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:17: ( phylonet_command )*
                loop5:
                do {
                    int alt5=2;
                    int LA5_0 = input.LA(1);

                    if ( ((LA5_0 >= BEGIN && LA5_0 <= DATA)||LA5_0==DIMENSIONS||(LA5_0 >= FORMAT && LA5_0 <= NCHAR)||(LA5_0 >= NETWORK && LA5_0 <= QUOTE)||(LA5_0 >= SYMBOLS && LA5_0 <= UTREE)||LA5_0==33||(LA5_0 >= 37 && LA5_0 <= 38)) ) {
                        alt5=1;
                    }


                    switch (alt5) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:106:17: phylonet_command
                        {
                            pushFollow(FOLLOW_phylonet_command_in_phylonet_block_body205);
                            phylonet_command();

                            state._fsp--;


                        }
                        break;

                        default :
                            break loop5;
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



    // $ANTLR start "data_block_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:108:1: data_block_body : DATA ';' (~ DIMENSIONS )* DIMENSIONS NTAX '=' identifier NCHAR '=' identifier ';' (~ FORMAT )* FORMAT DATATYPE '=' identifier SYMBOLS '=' QUOTE MISSING '=' ID GAP '=' ID ';' (~ MATRIX )* MATRIX ( identifier identifier )* ';' ;
    public final void data_block_body() throws RecognitionException {
        int numPairs = 0;
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:109:2: ( DATA ';' (~ DIMENSIONS )* DIMENSIONS NTAX '=' identifier NCHAR '=' identifier ';' (~ FORMAT )* FORMAT DATATYPE '=' identifier SYMBOLS '=' QUOTE MISSING '=' ID GAP '=' ID ';' (~ MATRIX )* MATRIX ( identifier identifier )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:109:4: DATA ';' (~ DIMENSIONS )* DIMENSIONS NTAX '=' identifier NCHAR '=' identifier ';' (~ FORMAT )* FORMAT DATATYPE '=' identifier SYMBOLS '=' QUOTE MISSING '=' ID GAP '=' ID ';' (~ MATRIX )* MATRIX ( identifier identifier )* ';'
            {
                match(input,DATA,FOLLOW_DATA_in_data_block_body222);

                match(input,37,FOLLOW_37_in_data_block_body224);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:110:10: (~ DIMENSIONS )*
                loop6:
                do {
                    int alt6=2;
                    int LA6_0 = input.LA(1);

                    if ( ((LA6_0 >= BEGIN && LA6_0 <= DEFAULT_INDICATOR)||(LA6_0 >= ELSE && LA6_0 <= 40)) ) {
                        alt6=1;
                    }


                    switch (alt6) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
                        {
                            if ( (input.LA(1) >= BEGIN && input.LA(1) <= DEFAULT_INDICATOR)||(input.LA(1) >= ELSE && input.LA(1) <= 40) ) {
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
                            break loop6;
                    }
                } while (true);


                match(input,DIMENSIONS,FOLLOW_DIMENSIONS_in_data_block_body242);

                match(input,NTAX,FOLLOW_NTAX_in_data_block_body244);

                match(input,39,FOLLOW_39_in_data_block_body246);

                pushFollow(FOLLOW_identifier_in_data_block_body248);
                identifier();

                state._fsp--;


                match(input,NCHAR,FOLLOW_NCHAR_in_data_block_body250);

                match(input,39,FOLLOW_39_in_data_block_body252);

                pushFollow(FOLLOW_identifier_in_data_block_body254);
                identifier();

                state._fsp--;


                match(input,37,FOLLOW_37_in_data_block_body256);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:112:3: (~ FORMAT )*
                loop7:
                do {
                    int alt7=2;
                    int LA7_0 = input.LA(1);

                    if ( ((LA7_0 >= BEGIN && LA7_0 <= END)||(LA7_0 >= GAP && LA7_0 <= 40)) ) {
                        alt7=1;
                    }


                    switch (alt7) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
                        {
                            if ( (input.LA(1) >= BEGIN && input.LA(1) <= END)||(input.LA(1) >= GAP && input.LA(1) <= 40) ) {
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
                            break loop7;
                    }
                } while (true);


                match(input,FORMAT,FOLLOW_FORMAT_in_data_block_body266);

                match(input,DATATYPE,FOLLOW_DATATYPE_in_data_block_body268);

                match(input,39,FOLLOW_39_in_data_block_body269);

                pushFollow(FOLLOW_identifier_in_data_block_body270);
                identifier();

                state._fsp--;


                match(input,SYMBOLS,FOLLOW_SYMBOLS_in_data_block_body272);

                match(input,39,FOLLOW_39_in_data_block_body273);

                match(input,QUOTE,FOLLOW_QUOTE_in_data_block_body274);

                match(input,MISSING,FOLLOW_MISSING_in_data_block_body276);

                match(input,39,FOLLOW_39_in_data_block_body277);

                match(input,ID,FOLLOW_ID_in_data_block_body278);

                match(input,GAP,FOLLOW_GAP_in_data_block_body280);

                match(input,39,FOLLOW_39_in_data_block_body281);

                match(input,ID,FOLLOW_ID_in_data_block_body282);

                match(input,37,FOLLOW_37_in_data_block_body284);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:114:3: (~ MATRIX )*
                loop8:
                do {
                    int alt8=2;
                    int LA8_0 = input.LA(1);

                    if ( ((LA8_0 >= BEGIN && LA8_0 <= ID_SET)||(LA8_0 >= MISSING && LA8_0 <= 40)) ) {
                        alt8=1;
                    }


                    switch (alt8) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
                        {
                            if ( (input.LA(1) >= BEGIN && input.LA(1) <= ID_SET)||(input.LA(1) >= MISSING && input.LA(1) <= 40) ) {
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
                            break loop8;
                    }
                } while (true);


                match(input,MATRIX,FOLLOW_MATRIX_in_data_block_body294);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:116:3: ( identifier identifier )*
                loop9:
                do {
                    int alt9=2;
                    int LA9_0 = input.LA(1);

                    if ( ((LA9_0 >= BEGIN && LA9_0 <= DATA)||LA9_0==DIMENSIONS||(LA9_0 >= FORMAT && LA9_0 <= ID)||(LA9_0 >= MATRIX && LA9_0 <= NCHAR)||(LA9_0 >= NETWORK && LA9_0 <= PHYLONET)||LA9_0==SYMBOLS||(LA9_0 >= TRANSLATE && LA9_0 <= UTREE)) ) {
                        alt9=1;
                    }


                    switch (alt9) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:116:4: identifier identifier
                        {
                            pushFollow(FOLLOW_identifier_in_data_block_body299);
                            identifier();

                            state._fsp--;


                            pushFollow(FOLLOW_identifier_in_data_block_body301);
                            identifier();

                            state._fsp--;


                            numPairs++;

                        }
                        break;

                        default :
                            break loop9;
                    }
                } while (true);


                stack.pushDataBlockBody(numPairs);

                match(input,37,FOLLOW_37_in_data_block_body310);

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
    // $ANTLR end "data_block_body"



    // $ANTLR start "phylonet_command"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:119:1: phylonet_command : ( identifier '=' )? ( phylonet_command_part )* ';' ;
    public final void phylonet_command() throws RecognitionException {
        boolean assignment = false;
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:2: ( ( identifier '=' )? ( phylonet_command_part )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:5: ( identifier '=' )? ( phylonet_command_part )* ';'
            {
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:5: ( identifier '=' )?
                int alt10=2;
                int LA10_0 = input.LA(1);

                if ( ((LA10_0 >= BEGIN && LA10_0 <= DATA)||LA10_0==DIMENSIONS||(LA10_0 >= FORMAT && LA10_0 <= ID)||(LA10_0 >= MATRIX && LA10_0 <= NCHAR)||(LA10_0 >= NETWORK && LA10_0 <= PHYLONET)||LA10_0==SYMBOLS||(LA10_0 >= TRANSLATE && LA10_0 <= UTREE)) ) {
                    int LA10_1 = input.LA(2);

                    if ( (LA10_1==39) ) {
                        alt10=1;
                    }
                }
                switch (alt10) {
                    case 1 :
                        // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:6: identifier '='
                    {
                        pushFollow(FOLLOW_identifier_in_phylonet_command336);
                        identifier();

                        state._fsp--;


                        match(input,39,FOLLOW_39_in_phylonet_command338);

                        assignment = true;

                    }
                    break;

                }


                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:45: ( phylonet_command_part )*
                loop11:
                do {
                    int alt11=2;
                    int LA11_0 = input.LA(1);

                    if ( ((LA11_0 >= BEGIN && LA11_0 <= DATA)||LA11_0==DIMENSIONS||(LA11_0 >= FORMAT && LA11_0 <= NCHAR)||(LA11_0 >= NETWORK && LA11_0 <= QUOTE)||(LA11_0 >= SYMBOLS && LA11_0 <= UTREE)||LA11_0==33||LA11_0==38) ) {
                        alt11=1;
                    }


                    switch (alt11) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:120:45: phylonet_command_part
                        {
                            pushFollow(FOLLOW_phylonet_command_part_in_phylonet_command345);
                            phylonet_command_part();

                            state._fsp--;


                        }
                        break;

                        default :
                            break loop11;
                    }
                } while (true);


                match(input,37,FOLLOW_37_in_phylonet_command348);

                stack.pushPhylonetCommand(assignment);

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:122:1: phylonet_command_part : ( identifier | ident_list |p= QUOTE |p= TAXON_SET_LIST |p= ID_SET | taxa_map );
    public final void phylonet_command_part() throws RecognitionException {
        Token p=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:123:2: ( identifier | ident_list |p= QUOTE |p= TAXON_SET_LIST |p= ID_SET | taxa_map )
            int alt12=6;
            switch ( input.LA(1) ) {
                case BEGIN:
                case DATA:
                case DIMENSIONS:
                case FORMAT:
                case GAP:
                case ID:
                case MATRIX:
                case MISSING:
                case NCHAR:
                case NETWORK:
                case NETWORKS:
                case NTAX:
                case PHYLONET:
                case SYMBOLS:
                case TRANSLATE:
                case TREE:
                case TREES:
                case UTREE:
                {
                    alt12=1;
                }
                break;
                case 33:
                {
                    alt12=2;
                }
                break;
                case QUOTE:
                {
                    alt12=3;
                }
                break;
                case TAXON_SET_LIST:
                {
                    alt12=4;
                }
                break;
                case ID_SET:
                {
                    alt12=5;
                }
                break;
                case 38:
                {
                    alt12=6;
                }
                break;
                default:
                    NoViableAltException nvae = new NoViableAltException("", 12, 0, input);
                    throw nvae;

            }

            switch (alt12) {
                case 1 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:123:4: identifier
                {
                    pushFollow(FOLLOW_identifier_in_phylonet_command_part360);
                    identifier();

                    state._fsp--;


                    stack.pushPhylonetCommandPartIdent();

                }
                break;
                case 2 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:124:4: ident_list
                {
                    pushFollow(FOLLOW_ident_list_in_phylonet_command_part379);
                    ident_list();

                    state._fsp--;


                    stack.pushPhylonetCommandPartIdentList();

                }
                break;
                case 3 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:125:4: p= QUOTE
                {
                    p=(Token)match(input,QUOTE,FOLLOW_QUOTE_in_phylonet_command_part395);

                    stack.pushPhylonetCommandPartQuote(p);

                }
                break;
                case 4 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:126:10: p= TAXON_SET_LIST
                {
                    p=(Token)match(input,TAXON_SET_LIST,FOLLOW_TAXON_SET_LIST_in_phylonet_command_part425);

                    stack.pushPhylonetCommandPartSetList(p);

                }
                break;
                case 5 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:127:4: p= ID_SET
                {
                    p=(Token)match(input,ID_SET,FOLLOW_ID_SET_in_phylonet_command_part440);

                    stack.pushPhylonetCommandPartIdSet(p);

                }
                break;
                case 6 :
                    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:128:4: taxa_map
                {
                    pushFollow(FOLLOW_taxa_map_in_phylonet_command_part454);
                    taxa_map();

                    state._fsp--;


                    stack.pushPhylonetCommandPartTaxaMap();

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:130:1: ident_list : s= '(' identifier ( ',' identifier )* ')' ;
    public final void ident_list() throws RecognitionException {
        Token s=null;

        int numIdentsInList = 1;
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:2: (s= '(' identifier ( ',' identifier )* ')' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:4: s= '(' identifier ( ',' identifier )* ')'
            {
                s=(Token)match(input,33,FOLLOW_33_in_ident_list486);

                pushFollow(FOLLOW_identifier_in_ident_list488);
                identifier();

                state._fsp--;


                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:21: ( ',' identifier )*
                loop13:
                do {
                    int alt13=2;
                    int LA13_0 = input.LA(1);

                    if ( (LA13_0==35) ) {
                        alt13=1;
                    }


                    switch (alt13) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:131:22: ',' identifier
                        {
                            match(input,35,FOLLOW_35_in_ident_list491);

                            pushFollow(FOLLOW_identifier_in_ident_list493);
                            identifier();

                            state._fsp--;


                            numIdentsInList++;

                        }
                        break;

                        default :
                            break loop13;
                    }
                } while (true);


                match(input,34,FOLLOW_34_in_ident_list500);

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:133:1: tree_assigment : tr= ( TREE | UTREE ) rich_newick_assignment ;
    public final void tree_assigment() throws RecognitionException {
        Token tr=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:134:2: (tr= ( TREE | UTREE ) rich_newick_assignment )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:134:4: tr= ( TREE | UTREE ) rich_newick_assignment
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


                pushFollow(FOLLOW_rich_newick_assignment_in_tree_assigment522);
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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:136:1: rich_newick_assignment : (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string ;
    public final void rich_newick_assignment() throws RecognitionException {
        Token d=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:137:2: ( (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:137:4: (d= DEFAULT_INDICATOR )? identifier '=' rich_newick_string
            {
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:137:5: (d= DEFAULT_INDICATOR )?
                int alt14=2;
                int LA14_0 = input.LA(1);

                if ( (LA14_0==DEFAULT_INDICATOR) ) {
                    alt14=1;
                }
                switch (alt14) {
                    case 1 :
                        // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:137:5: d= DEFAULT_INDICATOR
                    {
                        d=(Token)match(input,DEFAULT_INDICATOR,FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment536);

                    }
                    break;

                }


                pushFollow(FOLLOW_identifier_in_rich_newick_assignment539);
                identifier();

                state._fsp--;


                match(input,39,FOLLOW_39_in_rich_newick_assignment541);

                pushFollow(FOLLOW_rich_newick_string_in_rich_newick_assignment543);
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



    // $ANTLR start "skip_body"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:140:1: skip_body : ~ ( TREES | NETWORKS | PHYLONET | DATA ) ';' (~ END )* ;
    public final void skip_body() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:2: (~ ( TREES | NETWORKS | PHYLONET | DATA ) ';' (~ END )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:4: ~ ( TREES | NETWORKS | PHYLONET | DATA ) ';' (~ END )*
            {
                if ( input.LA(1)==BEGIN||(input.LA(1) >= DATATYPE && input.LA(1) <= NETWORK)||input.LA(1)==NTAX||(input.LA(1) >= QUOTE && input.LA(1) <= TREE)||(input.LA(1) >= UTREE && input.LA(1) <= 40) ) {
                    input.consume();
                    state.errorRecovery=false;
                }
                else {
                    MismatchedSetException mse = new MismatchedSetException(null,input);
                    throw mse;
                }


                match(input,37,FOLLOW_37_in_skip_body568);

                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:141:40: (~ END )*
                loop15:
                do {
                    int alt15=2;
                    int LA15_0 = input.LA(1);

                    if ( ((LA15_0 >= BEGIN && LA15_0 <= ELSE)||(LA15_0 >= FORMAT && LA15_0 <= 40)) ) {
                        alt15=1;
                    }


                    switch (alt15) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:
                        {
                            if ( (input.LA(1) >= BEGIN && input.LA(1) <= ELSE)||(input.LA(1) >= FORMAT && input.LA(1) <= 40) ) {
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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:144:1: rich_newick_string : (str=~ ( ';' ) )* ';' ;
    public final void rich_newick_string() throws RecognitionException {
        Token str=null;

        StringBuffer accum = new StringBuffer(); int line = -1; int col = -1;
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:145:2: ( (str=~ ( ';' ) )* ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:145:5: (str=~ ( ';' ) )* ';'
            {
                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:145:5: (str=~ ( ';' ) )*
                loop16:
                do {
                    int alt16=2;
                    int LA16_0 = input.LA(1);

                    if ( ((LA16_0 >= BEGIN && LA16_0 <= 36)||(LA16_0 >= 38 && LA16_0 <= 40)) ) {
                        alt16=1;
                    }


                    switch (alt16) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:145:6: str=~ ( ';' )
                        {
                            str=(Token)input.LT(1);

                            if ( (input.LA(1) >= BEGIN && input.LA(1) <= 36)||(input.LA(1) >= 38 && input.LA(1) <= 40) ) {
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


                match(input,37,FOLLOW_37_in_rich_newick_string614);

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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:150:1: identifier : s= ( TRANSLATE | DATA | DIMENSIONS | SYMBOLS | MISSING | MATRIX | GAP | TREE | UTREE | FORMAT | NETWORK | NCHAR | BEGIN | NETWORKS | NTAX | TREES | PHYLONET | ID ) ;
    public final void identifier() throws RecognitionException {
        Token s=null;

        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:2: (s= ( TRANSLATE | DATA | DIMENSIONS | SYMBOLS | MISSING | MATRIX | GAP | TREE | UTREE | FORMAT | NETWORK | NCHAR | BEGIN | NETWORKS | NTAX | TREES | PHYLONET | ID ) )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:151:4: s= ( TRANSLATE | DATA | DIMENSIONS | SYMBOLS | MISSING | MATRIX | GAP | TREE | UTREE | FORMAT | NETWORK | NCHAR | BEGIN | NETWORKS | NTAX | TREES | PHYLONET | ID )
            {
                s=(Token)input.LT(1);

                if ( (input.LA(1) >= BEGIN && input.LA(1) <= DATA)||input.LA(1)==DIMENSIONS||(input.LA(1) >= FORMAT && input.LA(1) <= ID)||(input.LA(1) >= MATRIX && input.LA(1) <= NCHAR)||(input.LA(1) >= NETWORK && input.LA(1) <= PHYLONET)||input.LA(1)==SYMBOLS||(input.LA(1) >= TRANSLATE && input.LA(1) <= UTREE) ) {
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
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:153:1: end_semi : END ';' ;
    public final void end_semi() throws RecognitionException {
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:153:9: ( END ';' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:153:11: END ';'
            {
                match(input,END,FOLLOW_END_in_end_semi709);

                match(input,37,FOLLOW_37_in_end_semi711);

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



    // $ANTLR start "taxa_map"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:1: taxa_map : s= '<' taxa_map_entry ( ';' taxa_map_entry )+ '>' ;
    public final void taxa_map() throws RecognitionException {
        Token s=null;

        int numKeys = 1;
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:155:36: (s= '<' taxa_map_entry ( ';' taxa_map_entry )+ '>' )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:2: s= '<' taxa_map_entry ( ';' taxa_map_entry )+ '>'
            {
                s=(Token)match(input,38,FOLLOW_38_in_taxa_map729);

                pushFollow(FOLLOW_taxa_map_entry_in_taxa_map730);
                taxa_map_entry();

                state._fsp--;


                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:22: ( ';' taxa_map_entry )+
                int cnt17=0;
                loop17:
                do {
                    int alt17=2;
                    int LA17_0 = input.LA(1);

                    if ( (LA17_0==37) ) {
                        alt17=1;
                    }


                    switch (alt17) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:157:23: ';' taxa_map_entry
                        {
                            match(input,37,FOLLOW_37_in_taxa_map733);

                            pushFollow(FOLLOW_taxa_map_entry_in_taxa_map735);
                            taxa_map_entry();

                            state._fsp--;


                            numKeys++;

                        }
                        break;

                        default :
                            if ( cnt17 >= 1 ) break loop17;
                            EarlyExitException eee =
                                    new EarlyExitException(17, input);
                            throw eee;
                    }
                    cnt17++;
                } while (true);


                match(input,40,FOLLOW_40_in_taxa_map741);

                stack.pushTaxaMap(numKeys, s);

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
    // $ANTLR end "taxa_map"



    // $ANTLR start "taxa_map_entry"
    // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:159:1: taxa_map_entry : identifier ':' identifier ( ',' identifier )* ;
    public final void taxa_map_entry() throws RecognitionException {
        int numValues = 1;
        try {
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:160:2: ( identifier ':' identifier ( ',' identifier )* )
            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:160:4: identifier ':' identifier ( ',' identifier )*
            {
                pushFollow(FOLLOW_identifier_in_taxa_map_entry757);
                identifier();

                state._fsp--;


                match(input,36,FOLLOW_36_in_taxa_map_entry759);

                pushFollow(FOLLOW_identifier_in_taxa_map_entry761);
                identifier();

                state._fsp--;


                // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:160:30: ( ',' identifier )*
                loop18:
                do {
                    int alt18=2;
                    int LA18_0 = input.LA(1);

                    if ( (LA18_0==35) ) {
                        alt18=1;
                    }


                    switch (alt18) {
                        case 1 :
                            // D:\\WorkDev\\Code\\Antlr\\Unstable\\PySON\\PySON.g:160:31: ',' identifier
                        {
                            match(input,35,FOLLOW_35_in_taxa_map_entry764);

                            pushFollow(FOLLOW_identifier_in_taxa_map_entry766);
                            identifier();

                            state._fsp--;


                            numValues++;

                        }
                        break;

                        default :
                            break loop18;
                    }
                } while (true);


                stack.pushTaxaMapEntry(numValues);

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
    // $ANTLR end "taxa_map_entry"

    // Delegated rules




    public static final BitSet FOLLOW_START_in_blocks52 = new BitSet(new long[]{0x0000000000000010L});
    public static final BitSet FOLLOW_block_in_blocks54 = new BitSet(new long[]{0x0000000000000010L});
    public static final BitSet FOLLOW_EOF_in_blocks59 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block67 = new BitSet(new long[]{0x0000000000100000L});
    public static final BitSet FOLLOW_networks_block_body_in_block69 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block71 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block76 = new BitSet(new long[]{0x0000000040000000L});
    public static final BitSet FOLLOW_trees_block_body_in_block78 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block83 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block88 = new BitSet(new long[]{0x0000000000400000L});
    public static final BitSet FOLLOW_phylonet_block_body_in_block90 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block92 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block97 = new BitSet(new long[]{0x0000000000000020L});
    public static final BitSet FOLLOW_data_block_body_in_block99 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block103 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_BEGIN_in_block108 = new BitSet(new long[]{0x000001FFBFAFFFD0L});
    public static final BitSet FOLLOW_skip_body_in_block110 = new BitSet(new long[]{0x0000000000000400L});
    public static final BitSet FOLLOW_end_semi_in_block122 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_NETWORKS_in_networks_block_body133 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_networks_block_body135 = new BitSet(new long[]{0x0000000000080002L});
    public static final BitSet FOLLOW_NETWORK_in_networks_block_body150 = new BitSet(new long[]{0x00000000F47BB9B0L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_networks_block_body152 = new BitSet(new long[]{0x0000000000080002L});
    public static final BitSet FOLLOW_TREES_in_trees_block_body172 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_trees_block_body174 = new BitSet(new long[]{0x00000000A0000002L});
    public static final BitSet FOLLOW_tree_assigment_in_trees_block_body189 = new BitSet(new long[]{0x00000000A0000002L});
    public static final BitSet FOLLOW_PHYLONET_in_phylonet_block_body201 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_phylonet_block_body203 = new BitSet(new long[]{0x00000062FCFBF932L});
    public static final BitSet FOLLOW_phylonet_command_in_phylonet_block_body205 = new BitSet(new long[]{0x00000062FCFBF932L});
    public static final BitSet FOLLOW_DATA_in_data_block_body222 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_data_block_body224 = new BitSet(new long[]{0x000001FFFFFFFFF0L});
    public static final BitSet FOLLOW_DIMENSIONS_in_data_block_body242 = new BitSet(new long[]{0x0000000000200000L});
    public static final BitSet FOLLOW_NTAX_in_data_block_body244 = new BitSet(new long[]{0x0000008000000000L});
    public static final BitSet FOLLOW_39_in_data_block_body246 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body248 = new BitSet(new long[]{0x0000000000020000L});
    public static final BitSet FOLLOW_NCHAR_in_data_block_body250 = new BitSet(new long[]{0x0000008000000000L});
    public static final BitSet FOLLOW_39_in_data_block_body252 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body254 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_data_block_body256 = new BitSet(new long[]{0x000001FFFFFFFFF0L});
    public static final BitSet FOLLOW_FORMAT_in_data_block_body266 = new BitSet(new long[]{0x0000000000000040L});
    public static final BitSet FOLLOW_DATATYPE_in_data_block_body268 = new BitSet(new long[]{0x0000008000000000L});
    public static final BitSet FOLLOW_39_in_data_block_body269 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body270 = new BitSet(new long[]{0x0000000004000000L});
    public static final BitSet FOLLOW_SYMBOLS_in_data_block_body272 = new BitSet(new long[]{0x0000008000000000L});
    public static final BitSet FOLLOW_39_in_data_block_body273 = new BitSet(new long[]{0x0000000000800000L});
    public static final BitSet FOLLOW_QUOTE_in_data_block_body274 = new BitSet(new long[]{0x0000000000010000L});
    public static final BitSet FOLLOW_MISSING_in_data_block_body276 = new BitSet(new long[]{0x0000008000000000L});
    public static final BitSet FOLLOW_39_in_data_block_body277 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_ID_in_data_block_body278 = new BitSet(new long[]{0x0000000000001000L});
    public static final BitSet FOLLOW_GAP_in_data_block_body280 = new BitSet(new long[]{0x0000008000000000L});
    public static final BitSet FOLLOW_39_in_data_block_body281 = new BitSet(new long[]{0x0000000000002000L});
    public static final BitSet FOLLOW_ID_in_data_block_body282 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_data_block_body284 = new BitSet(new long[]{0x000001FFFFFFFFF0L});
    public static final BitSet FOLLOW_MATRIX_in_data_block_body294 = new BitSet(new long[]{0x00000020F47BB930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body299 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_data_block_body301 = new BitSet(new long[]{0x00000020F47BB930L});
    public static final BitSet FOLLOW_37_in_data_block_body310 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_phylonet_command336 = new BitSet(new long[]{0x0000008000000000L});
    public static final BitSet FOLLOW_39_in_phylonet_command338 = new BitSet(new long[]{0x00000062FCFBF930L});
    public static final BitSet FOLLOW_phylonet_command_part_in_phylonet_command345 = new BitSet(new long[]{0x00000062FCFBF930L});
    public static final BitSet FOLLOW_37_in_phylonet_command348 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_phylonet_command_part360 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_ident_list_in_phylonet_command_part379 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_QUOTE_in_phylonet_command_part395 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_TAXON_SET_LIST_in_phylonet_command_part425 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_ID_SET_in_phylonet_command_part440 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_taxa_map_in_phylonet_command_part454 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_33_in_ident_list486 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_ident_list488 = new BitSet(new long[]{0x0000000C00000000L});
    public static final BitSet FOLLOW_35_in_ident_list491 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_ident_list493 = new BitSet(new long[]{0x0000000C00000000L});
    public static final BitSet FOLLOW_34_in_ident_list500 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_tree_assigment514 = new BitSet(new long[]{0x00000000F47BB9B0L});
    public static final BitSet FOLLOW_rich_newick_assignment_in_tree_assigment522 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_DEFAULT_INDICATOR_in_rich_newick_assignment536 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_rich_newick_assignment539 = new BitSet(new long[]{0x0000008000000000L});
    public static final BitSet FOLLOW_39_in_rich_newick_assignment541 = new BitSet(new long[]{0x000001FFFFFFFFF0L});
    public static final BitSet FOLLOW_rich_newick_string_in_rich_newick_assignment543 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_skip_body557 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_skip_body568 = new BitSet(new long[]{0x000001FFFFFFFBF2L});
    public static final BitSet FOLLOW_set_in_rich_newick_string593 = new BitSet(new long[]{0x000001FFFFFFFFF0L});
    public static final BitSet FOLLOW_37_in_rich_newick_string614 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_set_in_identifier628 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_END_in_end_semi709 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_end_semi711 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_38_in_taxa_map729 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_taxa_map_entry_in_taxa_map730 = new BitSet(new long[]{0x0000002000000000L});
    public static final BitSet FOLLOW_37_in_taxa_map733 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_taxa_map_entry_in_taxa_map735 = new BitSet(new long[]{0x0000012000000000L});
    public static final BitSet FOLLOW_40_in_taxa_map741 = new BitSet(new long[]{0x0000000000000002L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry757 = new BitSet(new long[]{0x0000001000000000L});
    public static final BitSet FOLLOW_36_in_taxa_map_entry759 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry761 = new BitSet(new long[]{0x0000000800000002L});
    public static final BitSet FOLLOW_35_in_taxa_map_entry764 = new BitSet(new long[]{0x00000000F47BB930L});
    public static final BitSet FOLLOW_identifier_in_taxa_map_entry766 = new BitSet(new long[]{0x0000000800000002L});

}