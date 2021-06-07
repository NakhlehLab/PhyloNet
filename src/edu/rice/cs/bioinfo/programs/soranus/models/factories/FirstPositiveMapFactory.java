package edu.rice.cs.bioinfo.programs.soranus.models.factories;

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
 * Time: 6:24 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FirstPositiveMapFactory<I,D>
{
    public Map<I,D> makeMapFromXMLFile(File xmlFile) throws ParserConfigurationException, IOException, SAXException {


        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

        // use the factory to take an instance of the document builder
        DocumentBuilder db = dbf.newDocumentBuilder();
        // parse using the builder to get the DOM mapping of the
        // XML file
        Document dom = db.parse(xmlFile);


        Element doc = dom.getDocumentElement();
        if(doc.getTagName().toLowerCase().equals("firstpositivedates"))
        {
            Map<I,D> resultMap = new HashMap<I, D>();

            NodeList firstPositiveElements = dom.getElementsByTagName("FirstPositive");
            for(int i = 0; i<firstPositiveElements.getLength(); i++)
            {
                Node element = firstPositiveElements.item(i);
                String dateText = element.getAttributes().getNamedItem("date").getNodeValue();
                String sourceIdText = element.getAttributes().getNamedItem("sourceId").getNodeValue();
                I id = makeId(sourceIdText);
                D date = makeDate(dateText);

                resultMap.put(id,date);
            }

            return resultMap;
        }

        throw new IllegalArgumentException();
    }

    protected abstract D makeDate(String dateText);

    protected abstract I makeId(String sourceIdText) ;
}
