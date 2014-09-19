#define FORCE_DEBUG

#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "base/WebServer.h"
#include <unistd.h>
#include <stdio.h>

#include <stdlib.h>

#include <unistd.h>
#include "util/SyncConn.h"
#include "base/StringUtil.h"
#include "visual/VisualAlignments.h"

const string sep = "%~!";

class ReqHandler : public IRequestHandler
{
public:
  ReqHandler() {
    //m_counter = 0;
    m_discard = 0;

    //m_input = "query.tangerine";
    //m_output = "tangerine.out";
    m_imageCounter = 0;
    m_resultCounter = 0;
    m_base = "./";
    m_pTrans = NULL;
    m_id_counter = 10;
  }

  void SetServer(const string &s) {
    m_server = s;
  }

  virtual void OnInit() {
    m_pWebServer->AddRouteThrough("temp.ps");
    m_pWebServer->AddRouteThrough("temp.pdf");
    m_pWebServer->AddRouteThrough("temp.jpeg");
    m_pWebServer->AddRouteThrough("temp.png");
    m_pWebServer->AddRouteThrough("semmel.jpg");
    m_pWebServer->AddRouteThrough("scilifelab.jpg");
    m_pWebServer->AddRouteThrough("uu.png");
    m_pWebServer->AddRouteThrough("bils_logo.png");
    m_pWebServer->AddRouteThrough("favicon.ico");

    m_pWebServer->AddRouteThrough("mojozoom.css");
    m_pWebServer->AddRouteThrough("mojozoom.js");
    m_pWebServer->AddRouteThrough("FAQ.html");
  }
  
  virtual int OnRequest(struct mg_connection *conn, const string & uri, const string & id) {
    cout << "OnRequest. URI: " << uri << " id: " << id << endl;
    cout << "ID: " << IDString(conn, id) << endl;
    if (uri == "/handle_post_request") {
      //m_counter++;
      //cout << "Counting: " << m_counter << endl;
      if (m_idsys.Counter(id) == 0 && strstr(id.c_str(), sep.c_str()) == NULL) {
	string data;
	m_pWebServer->GetVariable(conn, data, "input_1");
	cout << "Writing data: " << id << endl;
	WriteFile(data, conn, id);
	m_pWebServer->SendHTMLFile(conn, "please_wait.html");
	m_idsys.Inc(id);
      } else {
	cout << "Try SendReadResult w/ ID " << id << endl;
	if (ReadSendResult(conn, id)) {
	  //m_counter = 0;
	  m_idsys.Remove(id);
	} else {
	  m_pWebServer->SendHTMLFile(conn, "please_wait.html");
	  m_idsys.Inc(id);	 
	}
	    //m_pWebServer->SendHTMLFile(conn, "Results.html");     
      }
      return MG_TRUE;
    } else {
      // m_counter = 0;
      m_idsys.Add(id);
      m_pWebServer->SendHTMLFile(conn, "Tangerine_index.html", GetID());
      cout << "done" << endl;
      return MG_TRUE;
    }
  }
 
  virtual void OnAuth(struct mg_connection *conn) {

  }

  virtual void OnOther(struct mg_connection *conn) {

  }

  void ReadTaxaColors(const string & file) {
    
    if (file == "") {
      m_tax.Add("Metazoa;");
      m_tax.Add("Eukaryota;");
      m_tax.Add("Teleostei;");
      m_tax.Add("Archosauria;");
      m_tax.Add("Mammalia;");
      m_tax.Add("Cyanobacteria;");
      m_tax.Add("Amphibia;");
      m_tax.Add("Chordata;");
      m_tax.Add("Alveolata;");
      m_tax.Add("Viridiplantae;");
      m_tax.Add("Firmicutes;");
      m_tax.Add("Fungi;");
      m_tax.Add("Viruses;");
      m_tax.Add("Proteobacteria;");
      m_tax.Add("Embryophyta;");
      m_tax.Add("Actinobacteria;");
      m_tax.Add("Ascomycota;");
      m_tax.Add("Basidiomycota;");
      return;
    }

    FlatFileParser parser;
    parser.Open(file);

    while (parser.ParseLine()) {
      if (parser.GetItemCount() == 0)
	continue;
      m_tax.Add(parser.AsString(0));
    }

  }
  void ReadTaxonomy(const string & file) {
    m_tax.ReadTaxonomy(file);
  }

private:
  string GetID() {
    string s = Stringify(m_id_counter);
    m_id_counter++;
    return s;
  }

  void WriteFile(const string & data, struct mg_connection *conn, const string & id);
  bool ReadSendResult(struct mg_connection *conn, const string & id);
  int m_discard;  
  //int m_counter;
  IDSystem m_idsys;

  //string m_input;
  //string m_output;
  int m_imageCounter;
  int m_resultCounter;
  string m_base;
  TaxaColor m_tax;
  SCommTransmitter * m_pTrans;
  string m_server;
  int m_id_counter;
  
};

void ReqHandler::WriteFile(const string & data, struct mg_connection *conn, const string & id)
{
  m_discard = 0;

  /*
  FlatFileParser parser;
  parser.Open(m_output);
  while (parser.ParseLine()) {
    m_discard++;
  }
  cout << "Discarding lines: " << m_discard << endl;
  */


  //;SyncConnClient clie
  int n = strlen(data.c_str());
  if (n < 5)
    return;

  string input = IDString(conn, id) + ".input.fasta";
  string output = IDString(conn, id) + ".output";
  string done = output + ".done";
  string cmmd = "rm " + done;
  string server_message = input + " " + output;

  int rr = system(cmmd.c_str());

  FILE * pOut = fopen(input.c_str(), "w");
  //if (data[0] == ' ')
  //off = 1;
  // for (int i=0; i<strlen(buffer); i++) {
  //if (buffer[i] < 'A' || buffer[i] > 'Z')
  //  buffer[i] = ' ';
  //}

 
  bool b = false;
  for (int i=0; i<n; i++) {
    //cout << (int)data[i] << " " << data[i] << endl;
    //if (data[i] < 'A' || data[i] > 'Z')
    //continue;

 
   if (data[i] == 13 || data[i] == ' ')
      continue;

   if (!b && data[i] != '>')
     fprintf(pOut, ">UserQuery\n");
   b = true;

    fprintf(pOut, "%c", data[i]);
  }
  fprintf(pOut, "\n");
  fclose(pOut);

  
  m_pTrans = GetTransmitter();
  m_pTrans->SendWait(server_message.c_str());
  delete m_pTrans;
  m_pTrans = NULL;

}
 

string GetEnsID(const string & s)
{
  StringParser links;
  links.SetLine(s);
  if (links.GetItemCount() < 3 || links.AsString(0) != "Summary")
    return "";
  
  int i;
  StringParser ens;
  ens.SetLine(s, "|");
  if (ens.GetItemCount() <= 3)
    return "";
  return ens.AsString(3);
}

bool ReqHandler::ReadSendResult(struct mg_connection *conn, const string & idraw)
{
  string id = idraw;
  string filter;
  StringParser idparser;
  idparser.SetLine(id, sep);
  if (idparser.GetItemCount() == 2) {
    filter = idparser.AsString(1);
    id = idparser.AsString(0);
  } 

  string output = IDString(conn, id) + ".output";
  string done = output + ".done";
  FILE * p = fopen(done.c_str(), "r");
  if (p == NULL)
    return false;
  fclose(p);

 
  m_pWebServer->SendData(conn, "<script type=\"text/javascript\" src=\"mojozoom.js\"></script>");  
  m_pWebServer->SendData(conn, "<link type=\"text/css\" href=\"mojozoom.css\" rel=\"stylesheet\" />");  

  m_pWebServer->SendData(conn, "<html>  <head>");

  m_pWebServer->SendData(conn, "<meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\">");
  m_pWebServer->SendData(conn, "<title>Results</title>");
  m_pWebServer->SendData(conn, "</head>");
  m_pWebServer->SendData(conn, "<body>");
  m_pWebServer->SendData(conn, "<form method=\"POST\" action=\"/\">");

  m_pWebServer->SendData(conn, "<input type=\"submit\" value=\"Align more sequences\"/>");
  m_pWebServer->SendData(conn, "</form>");
 
  m_pWebServer->SendData(conn, "<br>Alignment Results: <pre>");
  int k = 0;

  // Sort & prepare the output
  string sort_cmmd = "../SortProtOutput ";
  if (filter != "") {
    sort_cmmd += " -filter ";
    sort_cmmd += "\"" + filter + "\"";    
  }
  sort_cmmd += " -i ";
  sort_cmmd += output + " > " + output + ".sorted";
  output += ".sorted";
  cout << "SORTING: " << sort_cmmd << endl;
  int rs = system(sort_cmmd.c_str());
    

  FlatFileParser parser;

  //string output = IDString(conn, id) + ".output";
  parser.Open(output);

  int i;
  svec<string> keep;

  string forImageFile = IDString(conn, id) + "forimage." /*+ Stringify(m_imageCounter)*/ + ".txt";
  string imageFile = IDString(conn, id) + "image." /*+ Stringify(m_imageCounter)*/ + ".temp.ps";
  string pdfFile = IDString(conn, id) + "image." /*+ Stringify(m_imageCounter)*/ + ".temp.pdf";
  //string jpegFile = IDString(conn, id) + "image." /*+ Stringify(m_imageCounter)*/ + ".temp.jpeg";
  string jpegFile = IDString(conn, id) + "image." /*+ Stringify(m_imageCounter)*/ + ".temp.png";
  m_imageCounter++;


  FILE * pForImage = fopen(forImageFile.c_str(), "w");
  while (parser.ParseLine()) {
    if (parser.GetItemCount() > 4 && parser.AsString(0) == "Best" && parser.AsString(1) == "hit") {
      cout << "Adding link: " << parser.AsString(3) << endl;
      string summary = parser.AsString(0) + " " + parser.AsString(1) + " " + parser.AsString(2);
    
      summary += "<a href=\"handle_post_request?";
      summary += id;
      summary += sep;
      summary += &(parser.AsString(3))[1];
      summary += "\"> ";

      summary += parser.AsString(3);
      summary += " </a>";
      for (i=4; i<parser.GetItemCount(); i++) {
	summary += " ";
	summary += parser.AsString(i);
      }

      summary += "<br>";


      m_pWebServer->SendData(conn, summary);
      if (filter == "")
	filter = &(parser.AsString(3))[1];
      continue;
    }

    //if (k < m_discard) {
    //k++;
    //continue;
    //}
    keep.push_back(parser.Line());
    fprintf(pForImage, "%s\n", parser.Line().c_str());
    //string send = parser.Line();
    //send += '\n';
    //m_pWebServer->SendData(conn, send);
    //cout << "Sending line: " << parser.Line() << endl;
  }
  fclose(pForImage);
  m_pWebServer->SendData(conn, "<br><br>");


  //string cmmd = "../runVisualAlignments -i " + forImageFile + " -o " + imageFile + " -j";
  //int ret = system(cmmd.c_str());
  PrintAlignment(m_tax, forImageFile, imageFile, true, "");

  // Make PDF in the background
  string pdf = "ps2pdf " + imageFile + /*" " + pdfFile +*/ " &";
  int rrr = system(pdf.c_str());


  string dataToSend = "<img alt=\"Visual results\" src=\"";
  dataToSend += jpegFile;
  dataToSend +=  "\" width=\"1000\" align=\"top\"";
  //---------------------------------------------------
  dataToSend += " data-zoomsrc=\"";
  dataToSend += jpegFile;
  dataToSend +=  "\" id=\"main_image\"";
    //---------------------------------------------------
  dataToSend += "><br>";
  m_pWebServer->SendData(conn, dataToSend);

  dataToSend = "<a href=\"";
  dataToSend += imageFile;
  dataToSend += "\">Postscript version of the image</a>  ";
  m_pWebServer->SendData(conn, dataToSend);

  dataToSend = "<a href=\"";
  dataToSend += pdfFile;
  dataToSend += "\">PDF version of the image</a><br><br>";
  m_pWebServer->SendData(conn, dataToSend);

  for (i=0; i<keep.isize(); i++) {
    string send = keep[i];
    send += '\n';
    m_pWebServer->SendData(conn, send);
    
    string ens = GetEnsID(send);
    if (ens != "") {
      dataToSend = "Search: <a target=\"_blank\" href=\"";
      dataToSend += "http://www.ensembl.org/Multi/Search/Results?q=";
      dataToSend += ens;
      dataToSend += ";y=0;site=ensembl_all;x=0\">";
      dataToSend += ens;
      dataToSend += " (Ensembl)</a> ";
      //http://www.ensembl.org/Multi/Search/Results?q=XP_003225938.1;y=0;site=ensembl_all;x=0
      m_pWebServer->SendData(conn, dataToSend);
      cout << dataToSend << endl;
      dataToSend = "<a target=\"_blank\" href=\"";
      dataToSend += "http://www.ncbi.nlm.nih.gov/gquery/?term=";
      dataToSend += ens;
      dataToSend += "\">";
      dataToSend += ens;
      dataToSend += " (NCBI)</a><br> ";
      m_pWebServer->SendData(conn, dataToSend);
      cout << dataToSend << endl;
    }
  }

  cout << "Lines: " << k << endl;
  m_pWebServer->SendData(conn, "/pre>");
  m_pWebServer->SendData(conn, "<br>");
  m_pWebServer->SendData(conn, "</body>");
  m_pWebServer->SendData(conn, "</html>");

  return true;
}

int main( int argc, char** argv )
{
  
  commandArg<string> fileCmmd("-t","taxonomy input file");
  commandArg<int> portCmmd("-port","listening port", 8080);
  commandArg<string> colCmmd("-col","taxa to color in", "");
  //commandArg<string> serverCmmd("-server","server name");
  commandLineParser P(argc,argv);
  P.SetDescription("Protein alignment web server.");
  P.registerArg(fileCmmd);
  P.registerArg(portCmmd);
  P.registerArg(colCmmd);
  //P.registerArg(serverCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string colName = P.GetStringValueFor(colCmmd);
  int port = P.GetIntValueFor(portCmmd);
  
  ReqHandler handler;
  handler.ReadTaxonomy(fileName);
  handler.ReadTaxaColors(colName);
  //handler.SetServer(server);

  WebServer server(&handler, port);

  
   
  return 0;
}
