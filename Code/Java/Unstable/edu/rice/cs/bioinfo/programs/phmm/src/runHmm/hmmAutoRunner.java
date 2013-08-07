package runHmm;

import java.io.File;
import java.io.FileInputStream;

public class hmmAutoRunner {

    /**
     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        FileInputStream fis = new FileInputStream(new File("ex/mouse/myautoinput.txt"));
        System.setIn(fis);

        runHmm.main(args);
    }

}
