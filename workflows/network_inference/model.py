import torch
import torch.nn as nn
import lightning as L

class ElementwiseLinear(nn.Module):
    def __init__(self, input_size: int) -> None:
        super(ElementwiseLinear, self).__init__()

        # w is the learnable weight of this layer module
        self.w = nn.Parameter(torch.rand(input_size), requires_grad=True)

    def forward(self, x: torch.tensor) -> torch.tensor:
        # simple elementwise multiplication
        return self.w * x

class EWlayer(nn.Module):
    def __init__(self, input_size):
        super(EWlayer, self).__init__()
        self.ew = ElementwiseLinear(input_size)

    def forward(self, x):
        x = self.ew(x)
        return x

class FClayer(nn.Module):
    def __init__(self, input_size, output_size):
        super(FClayer, self).__init__()
        self.fc = nn.Linear(input_size, output_size)

    def forward(self, x):
        x = self.fc(x)
        return x

## wraper for training with lightning
class LitWrapper(L.LightningModule):
    def __init__(self, model, criterion, learning_rate):
        super(LitWrapper, self).__init__()
        self.model = model
        self.criterion = criterion
        self.learning_rate = learning_rate
        self.losses = []

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_pred = self(x)
        loss = self.criterion(y_pred, y)
        self.log('train_loss', loss)
        self.losses.append(
            {"epoch": self.current_epoch, "batch": batch_idx, "train_loss": float(loss)}
        )
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        y_pred = self(x)
        loss = self.criterion(y_pred, y)
        self.log('val_loss', loss)
        self.losses.append(
            {"epoch": self.current_epoch, "batch": batch_idx, "val_loss": float(loss)}
        )
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)