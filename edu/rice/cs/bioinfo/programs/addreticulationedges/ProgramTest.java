package edu.rice.cs.bioinfo.programs.addreticulationedges;
import org.junit.*;

import java.math.BigDecimal;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/13/12
 * Time: 4:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProgramTest
{
    @Test
    public void testAddEdgesAndScale()
    {
        String result = new Program().addEdgesAndScale("(A:10,B:20)R;", 1, 1.0, new Random(13));
        Assert.assertEquals("((B:10)i2#H1:10::0.32885982758721066243623454283806495368480682373046875,(i2#H1:5::0.67114017241278933756376545716193504631519317626953125,A:5)i1:5)R;", result);

        result = new Program().addEdgesAndScale("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 1, 1.0, new Random(15));
        Assert.assertEquals("(((E:1)i2#H1:16::0.38062203043353715070651333007845096290111541748046875,(D:5,C:5)A:5)i1:5,(i2#H1:1::0.61937796956646284929348666992154903709888458251953125,F:2)B:20)R;", result);

        result = new Program().addEdgesAndScale("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 1, 2.0, new Random(15));
        Assert.assertEquals("(((E:2)i2#H1:32::0.38062203043353715070651333007845096290111541748046875,(D:10,C:10)A:10)i1:10,(i2#H1:2::0.61937796956646284929348666992154903709888458251953125,F:4)B:40)R;", result);

        result = new Program().addEdgesAndScale("((C:5,D:5)A:10,(E:2,F:2)B:20)R;", 2, 1.0, new Random(15));
        Assert.assertEquals("((((F:1)i4#H2:8::0.1689738364905790657388706677011214196681976318359375,(E:1)i2#H1:8::0.19175955892193741192386369220912456512451171875)i3:8,(D:5,C:5)A:5)i1:5,(i4#H2:1::0.8310261635094209342611293322988785803318023681640625,i2#H1:1::0.80824044107806258807613630779087543487548828125)B:20)R;", result);

    }

    @Test
    public void testSddEdgesAndEnforceHeight()
    {
           String result = new Program().addEdgesAndEnforceHeight("(A:10,B:10)R;", 1, 10, BigDecimal.ZERO, new Random(13));
           Assert.assertEquals("((B:2.5)i2#H1:7.5::0.32885982758721066243623454283806495368480682373046875,(i2#H1:2.5::0.67114017241278933756376545716193504631519317626953125,A:5)i1:5)R;", result);

           result = new Program().addEdgesAndEnforceHeight("(A:10,B:10)R;", 1, 20, BigDecimal.ZERO, new Random(13));
           Assert.assertEquals("((B:5.0)i2#H1:15.0::0.32885982758721066243623454283806495368480682373046875,(i2#H1:5.0::0.67114017241278933756376545716193504631519317626953125,A:10)i1:10)R;", result);

           result = new Program().addEdgesAndEnforceHeight("(A:10,B:10.1)R;", 1, 20, new BigDecimal(".1"), new Random(13));
           Assert.assertEquals("((B:5.0)i2#H1:15.0::0.32885982758721066243623454283806495368480682373046875,(i2#H1:5.0::0.67114017241278933756376545716193504631519317626953125,A:10)i1:10)R;", result);

           new Program().addEdgesAndEnforceHeight("(((10:4.014622972281253597849044467693602200597524642944335937500,9:4.014622972281253597849044467693602200597524642944335937500)i9:1.815308653021494893680767290788935497403144836425781250000,(8:4.181958544715561086195485529515281086787581443786621093750 (((7:0.993789958930273782999709197838456020690500736236572265625,6:0.993789958930273782999709197838456020690500736236572265625)i8:1.109331223292501170868593618479280848987400531768798828125,5:2.103121182222774953868302816317736869677901268005371093750)i7:0.276342699670375827241155519686799379996955394744873046875,4:2.379463881893150781109458336004536249674856662750244140625)i6:1.802494662822410305086027193510744837112724781036376953125)i5:1.647973148888034365718251450516618206165730953216552734375)i4:2.170068306396404069691961069565877551212906837463378906250,((3:0.492095239896266351730798049857185105793178081512451171875,2:0.492095239896266351730798049857185105793178081512451171875)i3:4.370760663544435479935508226390084018930792808532714843750,1:4.862855903440701831666306276247269124723970890045166015625)i2:3.137144096559297689939391773350507719442248344421386718750)i1;",
                                                            1, 8, new BigDecimal(".001"), new Random(13));
    }

    @Test
    public void existingHybridProbsUntouched()
    {
        String result = new Program().addEdgesAndScale("((C:1,X#1:1::.3)A:1,(D:1,(E:1)X#1:1::.7)B:1)R:1;", 1, 1.0, new Random(13));
        int i = 0;
    }

    @Test
    public void YunProblem()
    {
        String net = "(((5:1.0515605911113877104275074162192897673565390801161256370989414149898048822517893086114781908690929412841796875)i13#H1:5.8634052556904113117630696217852403652339946873372231108268297044400739481684325937749235890805721282958984375::0.9275174621507955574628567774198018014430999755859375,(((4:2.3794638818931513094565759243505589613972475159311969316836916089039945243488460846492671407759189605712890625,(i13#H1:1.0515605911113877104275074162192897673565390801161256370989414149898048822517893086114781908690929412841796875::0.0724825378492044425371432225801981985569000244140625,(6:0.9937899589302740036654080069741874491804795784712968510500843947346773232798256003661663271486759185791015625,7:0.9937899589302740036654080069741874491804795784712968510500843947346773232798256003661663271486759185791015625)i8:1.1093312232925014171896068254643920855325985817609544231477984352449324412237530168567900545895099639892578125)i7:0.2763426996703758886015610919119794266841693556989456574858087789243847598452674674263107590377330780029296875)i6:1.8024946628224107053202424793903591218107592710195501678970899731295658707797002762163174338638782501220703125,8:4.1819585447155620147768184037409180832080067869507470995807815820335603951285463608655845746397972106933593750)i5:1.6479731488880347316417982224200068848699475313730924787527405898208097967394536453866749070584774017333984375,(9:4.0146229722812544892744162708394734748839000490155283434933515139455412201741069111449178308248519897460937500,10:4.0146229722812544892744162708394734748839000490155283434933515139455412201741069111449178308248519897460937500)i9:1.8153086530214952967602599679375105347672405084881169574961471248689126589681563928024843335151672363281250000)i4:1.0850341531982022757719604118436051645125794491295091695922489475755086385522218961341422982513904571533203125)i11:1.0850341531982022757719604118436051645125794491295091695922489475755086385522218961341422982513904571533203125,(1:4.8628559034407029114372241630941022064296755383392689084877980143875077223558633932043449021875858306884765625,(2:0.4920952398962664609978911781121517178934187913130994902872931168083363295817633797923917882144451141357421875,3:0.4920952398962664609978911781121517178934187913130994902872931168083363295817633797923917882144451141357421875)i3:4.3707606635444364504393329849819504885362567470261694182005048975791713927741000134119531139731407165527343750)i2:3.1371440965592983865253132867540330906734376782435890090302220526178797466165804053161991760134696960449218750)i1;";

        net = new Program().addEdgesAndEnforceHeight(net, 1, 8, new BigDecimal(".001"), new Random(13));
        int i = 0;

    }
}
