# GNSS_RTK
Análises e códigos de desenvolvimento e montagem do GNSS RTK utilizando o Emlid Reach

* A pasta "Analises" contém os códigos em python para comparativo entre dados do GNSS utilizando RTK e do posicionamento absoluto por ponto simples.

* A pasta "NMEA_GNSS_to_GPS" contém um código para Arduino que converte a sentença GNGGA para GPGGA e corrige o checksum da mensagem. Isso possibilita a leitura dos dados do GNSS por qualquer dispositivo que só aceite sentenças GPGGA.
