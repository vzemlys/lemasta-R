freal1 <- c("y_r_sa","c_r_sa","x_r_sa-m_r_sa","i_r_sa" )
fnmr1 <- c("Bendrasis vidaus produktas (BVP)",
            "Namų ūkių vartojimo išlaidos",
            "Prekių ir paslaugų balansas",
            "Bendrojo pagrindinio kapitalo formavimas"
            )

tbreal1 <- produce.tb.sum(ftry$data,freal1,fnmr1,years=2005:2011)


