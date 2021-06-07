package edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser;

import edu.rice.cs.bioinfo.library.language.vaal._1_1.out.reading.parsers.VAALOutStringSplitParser;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.ParseException;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/31/13
 * Time: 12:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class DatafileRecogniser
{


    public KnownDatafileFormat recognise(File file)
    {
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        try {
            // use the factory to take an instance of the document builder
            DocumentBuilder db = dbf.newDocumentBuilder();
            // parse using the builder to get the DOM mapping of the
            // XML file
            Document dom = db.parse(file);


            Element doc = dom.getDocumentElement();
            if(doc.getTagName().toLowerCase().equals("sequencings"))
            {
                return new KnownDatafileFormat() {
                    public <R, T, E extends Exception> R execute(KnownDatafileFormatAlgo<R, T, E> algo, T input) throws E {
                        return algo.forSequencings(this, input);
                    }
                };
            }
            else if(doc.getTagName().toLowerCase().equals("traces"))
            {
                return new KnownDatafileFormat() {
                    public <R, T, E extends Exception> R execute(KnownDatafileFormatAlgo<R, T, E> algo, T input) throws E {
                        return algo.forTraces(this, input);
                    }
                };
            }
            else if(doc.getTagName().toLowerCase().equals("firstpositivedates"))
            {
                return new KnownDatafileFormat() {
                    public <R, T, E extends Exception> R execute(KnownDatafileFormatAlgo<R, T, E> algo, T input) throws E {
                        return algo.forFirstPositive(this, input);
                    }
                };
            }

        } catch (ParserConfigurationException pce) {
            int i = 0;
        } catch (SAXException se) {
            int i = 0;
        } catch (IOException ioe) {
            int i = 0;
        }

        try
        {
            new VAALOutStringSplitParser().parse(new FileInputStream(file));
            // if no exceptions
            return new KnownDatafileFormat()
            {
                public <R, T, E extends Exception> R execute(KnownDatafileFormatAlgo<R, T, E> algo, T input) throws E {
                    return algo.forVAALOut(this, input);
                }
            };
        }
        catch(FileNotFoundException e)
        {
        }
        catch(IOException e)
        {
        }
        catch(ParseException e)
        {
        }

        throw new IllegalArgumentException();
    }
}
