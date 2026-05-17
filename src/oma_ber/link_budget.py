"""Simple link-budget helpers for OMA at the photodiode input."""


def pd_input_oma_dbm(
    tx_oma_dbm: float,
    losses_db: list[float] | tuple[float, ...],
    margin_db: float = 0.0,
) -> float:
    """Calculate photodiode-input OMA in dBm from TX OMA and dB losses.

    tx_oma_dbm is the transmitter OMA in dBm. losses_db and margin_db are
    power penalties in dB.
    """
    if not isinstance(losses_db, list | tuple):
        msg = "losses_db must be a list or tuple of dB losses."
        raise ValueError(msg)
    if any(loss_db < 0 for loss_db in losses_db):
        msg = "losses_db values must be non-negative."
        raise ValueError(msg)
    if margin_db < 0:
        msg = "margin_db must be non-negative."
        raise ValueError(msg)

    return tx_oma_dbm - sum(losses_db) - margin_db
