class RedElectricaConPerdidasNegativas:
    def __init__(self, N_KM=400, N_molinos=150):
        # Kilómetros con pérdidas negativas
        self.kilometros = [KilometroConPerdidasNegativas(i, 10, 15, 0.85, 0.90, 1.5) 
                          for i in range(N_KM)]
        
        # Molinos (100 con Quijote, 50 sin Quijote)
        self.molinos_quijote = [MolinoConQuijote(i) for i in range(500)]
        self.molinos_sin_quijote = [MolinoSinQuijote(i) for i in range(250)]
        
        self.energia_total = 0.0
        self.potencia_total = 0.0
        
    def paso(self, t, viento):
        """Un paso de simulación con pérdidas negativas."""
        # 1. Molinos generan según viento
        P_molinos = 0
        for m in self.molinos_quijote:
            P_molinos += m.generar(viento) * 1.25  # Quijote +25%
        for m in self.molinos_sin_quijote:
            P_molinos += m.generar(viento) * 0.75  # Sin Quijote -25%
        
        # 2. Kilómetros descargan (con pérdidas negativas)
        P_KM = 0
        for km in self.kilometros:
            # Cada Kilómetro genera según su fase
            if km.stock_ALTA > 0:
                E_balance, _, _ = km.ciclo_con_patada(7)  # 7 fases = más potencia
                P_KM += E_balance / 0.1  # Potencia instantánea
            else:
                km.reset_potencial()
        
        # 3. Balance total
        P_total = P_molinos + P_KM
        self.energia_total += P_total * 0.1  # dt = 0.1s
        self.potencia_total = P_total
        
        return P_molinos, P_KM, P_total