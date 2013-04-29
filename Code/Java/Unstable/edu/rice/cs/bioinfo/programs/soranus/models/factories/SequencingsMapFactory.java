package edu.rice.cs.bioinfo.programs.soranus.models.factories;

import edu.rice.cs.bioinfo.programs.soranus.models.data.Sequencing;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/3/13
 * Time: 5:32 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SequencingsMapFactory<I>
{
    public Map<Sequencing,I> makeMapFromXMLFile(File xmlFile) throws ParserConfigurationException, IOException, SAXException {


        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

        // use the factory to take an instance of the document builder
        DocumentBuilder db = dbf.newDocumentBuilder();
        // parse using the builder to get the DOM mapping of the
        // XML file
        Document dom = db.parse(xmlFile);


        Element doc = dom.getDocumentElement();
        if(doc.getTagName().toLowerCase().equals("sequencings"))
        {
            Map<Sequencing,I> resultMap = new HashMap<Sequencing, I>();

            NodeList sequencingElements = dom.getElementsByTagName("Sequencing");
            for(int i = 0; i<sequencingElements.getLength(); i++)
            {
                Node element = sequencingElements.item(i);
                String sequenceText = element.getAttributes().getNamedItem("sequence").getNodeValue();
                String sourceIdText = element.getAttributes().getNamedItem("sourceId").getNodeValue();
                final I id = makeId(sourceIdText);

                resultMap.put(new Sequencing(sequenceText)
                {
                    @Override
                    public String toString()
                    {
                        return id.toString();
                    }
                }, id);
            }

            return resultMap;
        }

        throw new IllegalArgumentException();
    }

    protected abstract I makeId(String sourceIdText) ;
}
