package edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.parsers;

import edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.ir.LineDifferenceRecord;
import edu.rice.cs.bioinfo.library.language.vaal.out.reading.OutReadResult;
import edu.rice.cs.bioinfo.library.language.vaal.out.reading.csa.CSAError;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 3/21/13
 * Time: 1:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class VAALOutStringSplitParser
{
    public OutReadResult<List<LineDifferenceRecord>> parse(InputStream instream) throws IOException, ParseException
    {
        InputStreamReader reader = new InputStreamReader(instream);
        BufferedReader buffReader = new BufferedReader(reader);

        final LinkedList<LineDifferenceRecord> diffRecords = new LinkedList<LineDifferenceRecord>();

        try
        {
            String headerLine1 = buffReader.readLine();
            String headerLine2 = buffReader.readLine();
            String line = buffReader.readLine();

            if(headerLine1 == null || headerLine2 == null)
            {
                throw new ParseException("", -1);
            }

            while(line != null)
            {
                String[] lineElements = line.split("\\s+");

                if(lineElements.length < 2)
                    throw new ParseException("", -1);

                int contig = Integer.parseInt(lineElements[0]);
                int position = Integer.parseInt((lineElements[1]));

                String sample = null, ref = null;
                for(int i = 2; i<lineElements.length; i++)
                {
                    String[] keyAndValue = lineElements[i].split("=");
                    String key = keyAndValue[0];
                    String value = keyAndValue[1];

                    if(key.equals("sample"))
                    {
                        sample = value;
                    }
                    else if(key.equals("ref"))
                    {
                        ref = value;
                    }
                }

                if(sample == null || ref == null)
                    throw new ParseException("", -1);

                LineDifferenceRecord record = new LineDifferenceRecord(contig, position, sample, ref);
                diffRecords.add(record);

                line = buffReader.readLine();
            }

            return new OutReadResult<List<LineDifferenceRecord>>()
            {
                public List<LineDifferenceRecord> getDifferences() {
                    return diffRecords;
                }

                public Iterable<CSAError> getContextErrors() {
                    return null;
                }
            };

        }
        finally
        {
            buffReader.close();
            reader.close();
        }


    }
}
