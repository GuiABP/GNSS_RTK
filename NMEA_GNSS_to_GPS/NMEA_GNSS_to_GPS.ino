#include <SoftwareSerial.h>
char GPSBuffer[96]="";
char MsgNova[96]="";
char GPSIndex=0;
 
void setup() {
  // Tempo para inicialização do Emlid Reach e do ADCP
  delay(120000);
  // baud rate do ADCP é 38400
  Serial.begin(38400); 
}

void loop()
{
  // Checa a mensagem para converter
  CheckGPS();
}
 
void CheckGPS()
{
  int inByte;
  
  while (Serial.available() > 0)
  {
    inByte = Serial.read();

    // Inicio da mensagem
    if (inByte =='$')
    {
      GPSIndex = 0;
    }

    // Armazena todos os caracteres antes de chegar ao fim da mensagem
    if (inByte != '\r')
    {
      GPSBuffer[GPSIndex++] = inByte;
    }

    // Final da mensagem - Zera o index e processa a mensagem para conversao
    if (inByte == '\n')
    {
      ProcessGPSLine();
      GPSIndex = 0;
    }
  }
}

void ProcessGPSLine()
{
  int i = 0;
  
  // A seguir, é necessária realizar a correção do protocolo GNGGA e convertê-lo para GPGGA. O mesmo é válido para o GNVTG. Isso realiza
  // a conversão do protocolo para o ADCP aceitar o padrão GNSS disponibilizado pelo EMLID REACH.  

  //if verifica se a mensagem é do tipo desejado. Nesse caso é o GNVTG do GPS Emlid Reach.
  if ((GPSBuffer[1] == 'G') && (GPSBuffer[2] == 'N') && (GPSBuffer[3] == 'V') && (GPSBuffer[4] == 'T') && (GPSBuffer[5] == 'G'))
  {
    ConverteProtocolo();
  }
  
  //if verifica se a mensagem é do tipo desejado. Nesse caso é o GNGGA do GPS Emlid Reach.
  if ((GPSBuffer[1] == 'G') && (GPSBuffer[2] == 'N') && (GPSBuffer[3] == 'G') && (GPSBuffer[4] == 'G') && (GPSBuffer[5] == 'A'))
  {
    ConverteProtocolo();
  }
}

void ConverteProtocolo(){
  //Serial.println(GPSBuffer); //Print da mensagem original para conferir (deve ficar comentado no funcionamento real)
  
  GPSBuffer[2] = 'P';
  
  int cont = 1; //Contador para varrer a mensagem original
  int cont1 = 0; //Contador para incrementar a posição do vetor da mensagem convertida
  byte ChkSum; //Variável para armazenar o checksum hexadecimal calculado para corrigir o original
  
  while(GPSBuffer[cont]!='*'){ //Loop para eliminar o checksum original e deixar a mensagem sem ele
  //if(MsgOriginal[cont]!=','){
    MsgNova[cont1] = GPSBuffer[cont];
    cont1++;
  //}
  cont++;
  }

  ChkSum = checksum(MsgNova); //Calcula o checksum correto da mensagem
  //Serial.println(ChkSum, HEX); //Print para verificar o checksum calculado (deve ficar comentado no funcionamento real)

  //Prints a seguir recriam a mensagem de acordo com o protocolo original
  Serial.print('$');
  Serial.print(MsgNova);
  Serial.print('*');
  Serial.println(ChkSum, HEX);

  // Limpa os vetores para execução seguinte
  int j;
  for(j=0;j<96;j++){
    MsgNova[j]=0;
  }
  for(j=0;j<96;j++){
    GPSBuffer[j]=0;
  }
}

//Rotina a seguir usa a lógica XOR, de acordo com protocolo NMEA, para verificar o checksum certo
int checksum(const char *s) 
{
  int c = 0;
  while(*s)
      c ^= *s++;
  return c;
}
