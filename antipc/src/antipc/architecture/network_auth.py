"""Red + permisos — port conceptual de comunicacion hashcode (HMAC-SHA256)."""

from __future__ import annotations

import hashlib
import hmac
import secrets
from dataclasses import dataclass, field


@dataclass
class NetworkAuthGate:
    """
    Challenge-response HMAC (HASHTOOLCODE comunicacion hashcode).
    El runtime AntiPC no usa la red sin permiso verificado.
    """

    master_key: bytes = field(default_factory=lambda: secrets.token_bytes(32))
    _challenges: dict[str, bytes] = field(default_factory=dict)

    def derive_token(self, user_id: str) -> bytes:
        return hmac.new(self.master_key, user_id.encode(), hashlib.sha256).digest()

    def issue_challenge(self, user_id: str) -> tuple[str, bytes]:
        nonce = secrets.token_bytes(16)
        self._challenges[user_id] = nonce
        return user_id, nonce

    def verify_response(self, user_id: str, response: bytes) -> bool:
        nonce = self._challenges.pop(user_id, None)
        if nonce is None:
            return False
        expected = hmac.new(self.derive_token(user_id), nonce, hashlib.sha256).digest()
        return hmac.compare_digest(expected, response)

    def client_response(self, user_id: str, nonce: bytes) -> bytes:
        return hmac.new(self.derive_token(user_id), nonce, hashlib.sha256).digest()