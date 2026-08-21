# El ciclo de batería con pérdidas negativas
class KilometroConPerdidasNegativas:
    def __init__(self, id_modulo, m_peso, delta_h, eta_gen, eta_lift, E_perno):
        self.id = id_modulo
        self.m_peso = m_peso
        self.delta_h = delta_h
        self.eta_gen = eta_gen
        self.eta_lift = eta_lift
        self.E_perno = E_perno
        self.stock_ALTA = 10
        self.stock_BAJA = 0
        self.energia_generada = 0.0
        self.energia_consumida = 0.0
        self.ciclos = 0
        self.hurto_distancia = 0.0
        
    def ciclo_con_patada(self, fases_patada=5):
        """
        Ciclo completo con pérdidas negativas.
        fases_patada: 3, 5, 7 (más fases = más hurto)
        """
        # Factor de patada (más fases = más hurto)
        factor_patada = {3: 1.0, 5: 1.5, 7: 2.0}[fases_patada]
        
        # --- FASE 1: Bajada (Generación) ---
        if self.stock_ALTA > 0:
            # Bajada en recta: generas energía
            E_gen = self.eta_gen * self.m_peso * 9.81 * self.delta_h
            self.energia_generada += E_gen
            
            # --- FASE 2: Hurto de distancia ---
            # Robas distancia de recorrido
            self.hurto_distancia = self.delta_h * (1 - 1/factor_patada)
            E_hurto = self.m_peso * 9.81 * self.hurto_distancia
            self.energia_generada += E_hurto  # ¡Pérdida negativa!
            
            # --- FASE 3: Subida por flotación (casi gratis) ---
            # No consumes energía (flotabilidad)
            
            # --- FASE 4: Reset por perneo ---
            self.stock_ALTA -= 1
            self.stock_BAJA += 1
            self.energia_consumida += 4 * self.E_perno
            self.ciclos += 1
            
            # Balance del ciclo
            E_balance = self.energia_generada - self.energia_consumida
            return E_balance, self.energia_generada, self.energia_consumida
        else:
            # Modo pausa: resetear potencial
            if self.stock_BAJA > 0:
                E_reset = (self.m_peso * 9.81 * self.delta_h) / self.eta_lift
                self.energia_consumida += E_reset
                self.stock_BAJA -= 1
                self.stock_ALTA += 1
            return 0, 0, 0