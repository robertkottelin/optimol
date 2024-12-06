DEVELOPMENT = True


class Config:
    # Default CORS_ORIGINS for production
    CORS_ORIGINS = [
        "https://frontend1.com",
        "https://frontend2.com"
    ]

class DevelopmentConfig(Config):
    # CORS_ORIGINS for local development
    CORS_ORIGINS = ["http://localhost:3000"]

class ProductionConfig(Config):
    # CORS_ORIGINS for production with restricted origins
    CORS_ORIGINS = ["https://frontend1.com", "https://frontend2.com"]

class AllowAllConfig(Config):
    # "Allow all" option for CORS
    CORS_ORIGINS = "*"
