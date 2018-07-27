package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.ISMB2018;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.R;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPPseudoLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 1/27/18
 * Time: 5:12 PM
 * To change this template use File | Settings | File Templates.
 */

public class FigureRunningTime {
    private static void initNetHeights(Network<NetNodeInfo> network, double popSize) {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(network)) {
            if(node.isLeaf()) {
                node.setData(new NetNodeInfo(Utils.DEFAULT_NET_LEAF_HEIGHT));
                continue;
            }
            double height = Double.MAX_VALUE;
            for(NetNode<NetNodeInfo> child : node.getChildren()) {
                height = Math.min(height, child.getParentDistance(node) + child.getData().getHeight());
            }
            node.setData(new NetNodeInfo(height));
        }
        boolean setPopSize = true;
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(network)) {
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                node.setParentDistance(par, par.getData().getHeight() - node.getData().getHeight());
                if(node.getParentSupport(par) == node.NO_POP_SIZE && setPopSize) {
                    node.setParentSupport(par, popSize);
                }
            }
        }
        network.getRoot().setRootPopSize(popSize);
    }

    private static void adopt(NetNode<NetNodeInfo> par, NetNode<NetNodeInfo> child, double[] params) {
        par.adoptChild(child, par.getData().getHeight() - child.getData().getHeight());
        child.setParentProbability(par, params[0]);
        child.setParentSupport(par, params[1]);
    }

    private static double[] getParameters(NetNode<NetNodeInfo> par, NetNode<NetNodeInfo> child) {
        return new double[] {child.getParentProbability(par), child.getParentSupport(par)};
    }

    private static void addReticulation(NetNode<NetNodeInfo> v1, NetNode<NetNodeInfo> v2,
                                     NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                                     NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6,
                                     double gamma, double popSize) {

        double[] paramV3V4 = getParameters(v3, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        double a = Math.min(paramV3V4[1], paramV5V6[1]);
        double b = Math.max(paramV3V4[1], paramV5V6[1]);

        v3.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v1, new double[] {NetNode.NO_PROBABILITY, popSize} );
        adopt(v1, v4, paramV3V4);
        adopt(v5, v2, new double[] {1.0-gamma, popSize} );
        adopt(v2, v6, paramV5V6);

        adopt(v1, v2, new double[] {gamma, popSize} );
    }

    private static void addRandomReticulation(Network<NetNodeInfo> network) {
        NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5, _v6;
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = Networks.getAllEdges(network);
        int numEdges = edges.size();
        int numRetiNodes = network.getReticulationCount();

        Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge1, edge2;
        edge1 = edges.get(Randomizer.getRandomInt(numEdges));
        do {
            edge2 = edges.get(Randomizer.getRandomInt(numEdges));
        } while (edge1 == edge2);

        _v3 = edge1.Item2;
        _v4 = edge1.Item1;
        _v5 = edge2.Item2;
        _v6 = edge2.Item1;

        double t3 = _v3.getData().getHeight();
        double t4 = _v4.getData().getHeight();
        double l1 = t3 - t4;
        double t1 = t4 + 0.5 * l1;

        double t5 = _v5.getData().getHeight();
        double t6 = _v6.getData().getHeight();
        double l2 = t5 - t6;
        double t2 = t6 + 0.5 * l2;

        double gamma = 0.5;

        _v1 = new BniNetNode<>();
        _v1.setData(new NetNodeInfo(t1));
        _v2 = new BniNetNode<>();
        _v2.setData(new NetNodeInfo(t2));

        if(t1 > t2)
            addReticulation(_v1, _v2, _v3, _v4, _v5, _v6, gamma, network.getRoot().getRootPopSize()); // v3 v4 t1 v5 v6 t2
        else
            addReticulation(_v2, _v1, _v5, _v6, _v3, _v4, gamma, network.getRoot().getRootPopSize());
    }

    private static boolean isNetworkValid(Network<NetNodeInfo> network){
        Set<String> leaves = new HashSet<>();
        for(NetNode node : network.getLeaves()) {
            leaves.add(node.getName());
        }

        int count = 0;
        for(Object leaf: network.getLeaves()){
            if(leaves.contains(((NetNode)leaf).getName())){
                count++;
            }
            else{
                return false;
            }
        }
        if(count!=leaves.size())return false;
        if(!Networks.isDisconnectedNetwork(network,null))return false;
        for(Object node: Networks.postTraversal(network)){
            double totalProb = 0;
            for (Object parent : ((NetNode) node).getParents()) {
                totalProb += ((NetNode) node).getParentProbability((NetNode) parent);
            }
            if(((NetNode)node).getChildCount()==1 && ((NetNode)node).getParentCount()<2){
                return false;
            }
            if(totalProb!=NetNode.NO_PROBABILITY && ((NetNode)node).isNetworkNode()){
                if(Math.abs(totalProb - 1) > 0.00001) {
                    throw new RuntimeException(network.toString());
                }
            }
            else if(!((NetNode)node).isRoot()){
                if(totalProb != NetNode.NO_PROBABILITY){
                    throw new RuntimeException(network.toString());
                }
            }
        }
        return true;
    }

    public static void go() {
        Utils._NUM_THREADS = 8;

        int numSites = 10000;
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;

        List<Network<NetNodeInfo>> starters = new ArrayList<>();
        starters.add(Networks.readNetwork("(1:22.48683929,((2:8.05934715,(((3:1.87509346,4:1.87509346):1.42588234,((5:1.20218468,6:1.20218468):0.04821014,(7:0.07337952,8:0.07337952):1.17701530):2.05058098):0.49372292,9:3.79469872):4.26464844):0.01494217,10:8.07428932):14.41254997);"));
        starters.add(Networks.readNetwork("(((((1:0.60886765,2:0.60886765):0.21164322,3:0.82051086):0.46471024,4:1.28522110):1.79495811,5:3.08017921):25.19387436,((((6:0.43833160,7:0.43833160):5.23575974,8:5.67409134):0.11312294,(9:3.90265083,10:3.90265083):1.88456345):8.05934715,(((((11:0.33256912,12:0.33256912):2.67498207,13:3.00755119):0.01170921,14:3.01926041):4.01834869,((15:0.44890976,16:0.44890976):5.07913971,17:5.52804947):1.50955963):2.54430389,(18:2.57087708,(19:0.31581497,20:0.31581497):2.25506210):7.01103592):4.26464844):14.42749214);"));
        starters.add(Networks.readNetwork("((((1:9.57980728,((2:1.44622040,3:1.44622040):6.36743164,(4:6.85935211,5:6.85935211):0.95429993):1.76615524):0.67635345,((6:0.38061523,7:0.38061523):3.03581238,(8:2.60424805,(9:1.79208374,10:1.79208374):0.81216431):0.81217957):6.83973312):1.79495811,((((11:2.37528229,12:2.37528229):0.48159027,(((13:0.42972183,(14:0.41870880,15:0.41870880):0.01101303):0.22890091,16:0.65862274):1.65975571,17:2.31837845):0.53849411):3.32643890,((18:0.55014038,19:0.55014038):2.61812210,20:3.16826248):3.01504898):2.31235504,(((21:0.21119690,22:0.21119690):3.10795212,23:3.31914902):4.01662827,(24:6.64775085,(25:1.64096069,26:1.64096069):5.00679016):0.68802643):1.15988922):3.55545235):25.19387436,((((27:8.21447754,28:8.21447754):6.43055344,((29:2.92562866,30:2.92562866):0.27816391,31:3.20379257):11.44123840):0.11312294,32:14.75815392):8.05934715,((33:16.00854874,(34:2.73012161,(35:2.52329254,36:2.52329254):0.20682907):13.27842712):2.54430389,((((37:0.18454742,38:0.18454742):3.54642105,((39:1.20355988,40:1.20355988):0.30421448,41:1.50777435):2.22319412):5.00809479,(42:4.41995621,43:4.41995621):4.31910706):2.80275345,(44:9.28675461,((((45:0.44469833,46:0.44469833):0.01347351,47:0.45817184):3.35312271,(48:0.04359436,49:0.04359436):3.76770020):5.07624435,50:8.88753891):0.39921570):2.25506210):7.01103592):4.26464844):14.42749214);"));
        starters.add(Networks.readNetwork("((((1:15.27171326,((2:7.13812637,((3:2.62088776,(4:2.18691254,(((5:0.02330399,6:0.02330399):1.18935776,7:1.21266174):0.16494370,8:1.37760544):0.80930710):0.43397522):3.07101822,(((((9:0.01859665,10:0.01859665):2.29739761,11:2.31599426):1.05616760,12:3.37216187):1.35832977,13:4.73049164):0.66957474,14:5.40006638):0.29183960):1.44622040):6.36743164,(((15:0.89663315,16:0.89663315):3.41148376,(17:3.87561417,(((18:0.71570206,19:0.71570206):1.21697998,(20:0.73489761,21:0.73489761):1.19778442):0.99386978,(22:0.35624313,23:0.35624313):2.57030869):0.94906235):0.43250275):8.24314117,(24:2.09769440,25:2.09769440):10.45356369):0.95429993):1.76615524):0.67635345,(((26:4.73426819,27:4.73426819):1.33825302,28:6.07252121):3.03581238,((29:4.63269806,(30:1.24493027,31:1.24493027):3.38776779):3.66345596,(32:7.48398972,33:7.48398972):0.81216431):0.81217957):6.83973312):1.79495811,((((34:8.06718826,(35:2.99277115,(36:1.45819473,37:1.45819473):1.53457642):5.07441711):0.48159027,((38:6.35052872,(39:0.66706848,40:0.66706848):5.68346024):1.65975571,41:8.01028442):0.53849411):3.32643890,(42:8.86016846,43:8.86016846):3.01504898):2.31235504,((44:9.01105499,(45:4.64778519,(((46:0.16878510,47:0.16878510):1.51924515,(48:0.63956070,49:0.63956070):1.04846954):2.70622253,(50:2.32404709,51:2.32404709):2.07020569):0.25353241):4.36326981):4.01662827,(((52:1.80519104,53:1.80519104):2.36446381,((((54:0.27454376,55:0.27454376):0.60758209,(56:0.58839798,57:0.58839798):0.29372787):0.43958664,58:1.32171249):1.88119125,59:3.20290375):0.96675110):8.17000198,(60:0.39506912,61:0.39506912):11.94458771):0.68802643):1.15988922):3.55545235):25.19387436,((((62:13.90638351,(63:2.25793457,64:2.25793457):11.64844894):6.43055344,((65:0.71003723,66:0.71003723):0.45694733,67:1.16698456):19.16995239):0.11312294,68:20.45005989):8.05934715,((((69:1.93278885,70:1.93278885):1.39759445,(71:3.24524689,((72:0.03075409,73:0.03075409):0.66006470,74:0.69081879):2.55442810):0.08513641):18.37007141,((75:1.92652893,76:1.92652893):6.49549866,(77:0.21542358,78:0.21542358):8.20660400):13.27842712):2.54430389,((((79:5.87645340,80:5.87645340):3.54642105,((81:4.74730682,82:4.74730682):2.14815903,83:6.89546585):2.52740860):5.00809479,((84:3.70026398,85:3.70026398):6.41159821,((86:3.29062271,(87:0.20001984,88:0.20001984):3.09060287):2.00871277,(89:2.51729965,90:2.51729965):2.78203583):4.81252670):4.31910706):2.80275345,(((91:6.13660431,92:6.13660431):0.01347351,(93:1.34955597,(94:0.20747757,95:0.20747757):1.14207840):4.80052185):3.35312271,((96:4.73130798,(97:1.25988007,98:1.25988007):3.47142792):1.00419235,(99:1.21600342,100:1.21600342):4.51949692):3.76770020):7.73052216):7.01103592):4.26464844):14.42749214);"));
        starters.add(Networks.readNetwork("((((((1:0.10647202,2:0.10647202):13.96485138,((((3:0.76581573,4:0.76581573):1.33012009,5:2.09593582):7.45814896,((6:4.71556854,7:4.71556854):4.40454102,(8:0.57367325,9:0.57367325):8.54643631):0.43397522):3.07101822,(((((10:5.80704880,11:5.80704880):3.44214249,(12:2.08835983,13:2.08835983):7.16083145):1.05616760,(((14:0.73749542,(15:0.71877670,(16:0.37144470,17:0.37144470):0.34733200):0.01871872):1.21592331,(18:0.94195938,19:0.94195938):1.01145935):2.66542816,((20:0.07100677,21:0.07100677):3.26436234,(22:1.33276367,23:1.33276367):2.00260544):1.28347778):5.68651199):1.35832977,(((24:0.91224670,25:0.91224670):3.71193314,(26:3.26011658,27:3.26011658):1.36406326):1.17564774,(28:2.85452271,29:2.85452271):2.94530487):5.86386108):0.66957474,((30:0.41208649,(31:0.20080948,32:0.20080948):0.21127701):3.31270599,((33:0.19566345,34:0.19566345):2.77031708,35:2.96598053):0.75881195):8.60847092):0.29183960):1.44622040):6.36743164,((36:11.24131393,(((37:4.20320892,(38:3.36856079,39:3.36856079):0.83464813):0.82178497,40:5.02499390):5.78381729,((((((41:0.22285461,42:0.22285461):0.86232758,43:1.08518219):0.88026810,(44:0.33981705,45:0.33981705):1.62563324):5.68344879,(46:1.81651306,47:1.81651306):5.83238602):1.21697998,((48:3.05837631,49:3.05837631):4.60971832,((50:0.84201050,51:0.84201050):1.10247421,52:1.94448471):5.72360992):1.19778442):0.99386978,(((53:0.53849030,54:0.53849030):2.41473389,55:2.95322418):4.33621597,(56:0.05718994,57:0.05718994):7.23225021):2.57030869):0.94906235):0.43250275):8.24314117,(((58:3.65450287,59:3.65450287):2.96479034,(60:1.16104126,61:1.16104126):5.45825195):2.41159821,62:9.03089142):10.45356369):0.95429993):2.44250870,(63:16.04153061,(((64:4.57080078,65:4.57080078):6.99509430,(66:8.17812729,((67:0.47724533,68:0.47724533):0.61637497,69:1.09362030):7.08450699):3.38776779):3.66345596,70:15.22935104):0.81217957):6.83973312):1.79495811,((((71:15.00038528,(72:9.92596817,(73:0.32448196,74:0.32448196):9.60148621):5.07441711):0.48159027,((75:2.39135742,76:2.39135742):5.20890808,77:7.60026550):7.88171005):3.32643890,(78:6.09919739,(79:2.21403122,((80:0.56494141,81:0.56494141):1.30165482,(82:1.76055145,83:1.76055145):0.10604477):0.34743500):3.88516617):12.70921707):2.31235504,((((84:2.13806534,85:2.13806534):1.40977478,(86:1.58432388,87:1.58432388):1.96351624):12.39641190,(((88:0.60158157,(89:0.57497406,90:0.57497406):0.02660751):1.67653656,(91:1.35011292,92:1.35011292):0.92800522):9.30286407,(((93:5.65525818,94:5.65525818):1.44672394,95:7.10198212):1.51924515,96:8.62122726):2.95975494):4.36326981):4.01662827,((((97:3.21546555,(98:2.26839828,99:2.26839828):0.94706726):1.15248871,(100:1.40347672,101:1.40347672):2.96447754):6.73489761,((((102:0.14767838,103:0.14767838):7.66764450,(((104:1.38702393,(105:0.52935791,106:0.52935791):0.85766602):3.49468231,(107:1.13793945,108:1.13793945):3.74376678):2.63988876,(109:1.84632111,110:1.84632111):5.67527390):0.29372787):0.43958664,((111:0.65491486,112:0.65491486):0.96013260,113:1.61504745):6.63986206):1.88119125,(114:5.90146255,115:5.90146255):4.23463821):0.96675110):8.17000198,(((116:2.96414948,((117:0.07111359,118:0.07111359):1.40026855,119:1.47138214):1.49276733):2.31935501,(120:4.49320602,(121:0.89734268,122:0.89734268):3.59586334):0.79029846):0.30289078,123:5.58639526):13.68645859):0.68802643):1.15988922):3.55545235):25.19387436,((((((124:3.16495895,125:3.16495895):1.47830963,(126:3.34379196,127:3.34379196):1.29947662):16.19631195,(128:9.19113159,(129:6.13312912,130:6.13312912):3.05800247):11.64844894):6.43055344,((131:6.00009155,132:6.00009155):2.10009003,133:8.10018158):19.16995239):0.11312294,(134:1.68377686,(135:1.30065918,136:1.30065918):0.38311768):25.69948006):8.05934715,((((137:1.59856033,(138:0.81821442,139:0.81821442):0.78034592):8.66501999,(140:7.62401581,((((141:2.08995819,((142:0.21659088,143:0.21659088):0.09817505,144:0.31476593):1.77519226):0.13797760,145:2.22793579):0.26736069,(146:0.44220734,147:0.44220734):2.05308914):1.27450943,148:3.76980591):3.85420990):2.63956451):18.37007141,(149:15.35522461,(((150:0.17911148,151:0.17911148):6.71741867,(152:0.12313461,153:0.12313461):6.77339554):0.25209045,((154:0.07213974,155:0.07213974):3.36159897,((156:0.40941620,157:0.40941620):0.59955597,158:1.00897217):2.42476654):3.71488190):8.20660400):13.27842712):2.54430389,((((159:12.80965042,(160:0.41163254,161:0.41163254):12.39801788):3.54642105,((((162:2.86267853,163:2.86267853):3.19353485,((164:3.48155594,(165:2.77882004,(166:2.60213470,(167:1.24262619,168:1.24262619):1.35950851):0.17668533):0.70273590):0.76873779,169:4.25029373):1.80591965):5.62429047,(170:1.43576431,171:1.43576431):10.24473953):2.14815903,((((172:0.14590454,173:0.14590454):1.10869217,174:1.25459671):1.99287415,(175:2.18420029,(176:0.63695145,177:0.63695145):1.54724884):1.06327057):1.44983292,((178:4.20032883,179:4.20032883):0.20531845,180:4.40564728):0.29165649):9.13135910):2.52740860):5.00809479,((181:10.63346100,((182:2.58219528,183:2.58219528):2.67586899,((184:0.53560257,185:0.53560257):1.15195847,186:1.68756104):3.57050323):5.37539673):6.41159821,(187:1.82694626,188:1.82694626):15.21811295):4.31910706):2.80275345,((189:13.08327484,((190:2.01992798,191:2.01992798):0.46839523,(192:0.34364700,193:0.34364700):2.14467621):10.59495163):3.35312271,((((194:1.88154221,(195:0.32210159,196:0.32210159):1.55944061):2.10092163,197:3.98246384):7.68204117,198:11.66450500):1.00419235,(199:8.14920044,200:8.14920044):4.51949692):3.76770020):7.73052216):7.01103592):4.26464844):14.42749214);"));

        int numTaxa[] = {10, 20, 50, 100};
        int maxReti = 4;
        boolean usePseudoLikelihood = false;

        double timeConsumption[][] = new double[numTaxa.length][maxReti + 1];
        for(int i = 0 ; i < numTaxa.length ; i++) {
            Map<String, String> alleles2species = new HashMap<>();
            for(int j = 1 ; j <= numTaxa[i]; j++) {
                String s = String.format("%d", j);
                alleles2species.put(s, s);
            }

            Network<NetNodeInfo> currentNetwork = starters.get(i);
            System.out.println("# of taxa = " + numTaxa[i]);

            initNetHeights(currentNetwork, 0.01);

            for(int numReti = 0 ; numReti <= maxReti ; numReti++) {
                if(numReti > 0) {
                    addRandomReticulation(currentNetwork);
                    if(Networks.hasCycle(currentNetwork) || !isNetworkValid(currentNetwork)) {
                        throw new RuntimeException("invalid network: i=" + i + " numReti=" + numReti + " " + currentNetwork);
                    }
                }
                System.out.println("# of reticulations = " + currentNetwork.getReticulationCount());

                SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
                simulator._diploid = false;
                Map<String, String> onesnp = simulator.generateSNPs(currentNetwork, null, numSites, !useOnlyPolymorphic);

                List<Alignment> alns = new ArrayList<>();
                Alignment aln = new Alignment(onesnp);
                aln._diploid = false;
                alns.add(aln);

                if(usePseudoLikelihood) {
                    SNAPPPseudoLikelihood pseudoLikelihood = new SNAPPPseudoLikelihood(alleles2species, alns, alns.get(0)._diploid);

                    long start = System.currentTimeMillis();
                    double ll = pseudoLikelihood.computeSNAPPPseudoLogLikelihoodMT(currentNetwork, alleles2species, BAGTRModel);
                    long end = System.currentTimeMillis();
                    System.out.println("Time(s): " + (end - start) / 1000.0);
                    timeConsumption[i][numReti] = (end - start) / 1000.0;

                    System.out.println("likelihood = " + ll);
                } else {
                    SNAPPLikelihood.timeSavingMode = true;
                    aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns);
                    Network cloneNetwork = currentNetwork.clone();
                    R.maxLineages = aln._RPatterns.keySet().iterator().next().sumLineages();
                    long start = System.currentTimeMillis();
                    double ll = SNAPPLikelihood.computeSNAPPLikelihood(cloneNetwork, aln._RPatterns, BAGTRModel);
                    long end = System.currentTimeMillis();
                    System.out.println("Time(s): " + (end - start) / 1000.0);
                    timeConsumption[i][numReti] = (end - start) / 1000.0;
                    System.out.println("likelihood = " + ll);
                }
            }
        }

        for(int i = 0 ; i < numTaxa.length ; i++) {
            for(int numReti = 0 ; numReti <= maxReti ; numReti++) {
                System.out.print(timeConsumption[i][numReti] + " ");
            }
            System.out.println();
        }
    }

    public static void main(String[] args) {
        go();
    }
}

