# This container is a basic web frontend, we do the heavy-lifting in the
# main DiffDock image.
FROM rbgcsail/diffdock

USER $APPUSER
WORKDIR $HOME/$DIR_NAME

# Precompute series for SO(2) and SO(3) groups
RUN micromamba run -n ${ENV_NAME} python utils/precompute_series.py

# Expose port for web service
EXPOSE 7860

# Run app under micromamba environment
CMD ["sh", "-c", "micromamba run -n ${ENV_NAME} python3 app/main.py"]
